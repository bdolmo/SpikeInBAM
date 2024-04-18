#include "BamReader.h"
#include "BamRecord.h"
#include "BamWriter.h"
#include "RefFasta.h"
#include <iostream>
#include <fstream>
#include <htslib/sam.h>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include "sw.cpp"
#include <random>
#include "SnvCaller.cpp"
#include "Simulations.cpp"
#include <chrono>
#include <iomanip>
#include <libgen.h>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


std::string getBasename(const std::string& path) {
    size_t lastSlash = path.find_last_of("/\\");
    if (lastSlash != std::string::npos)
        return path.substr(lastSlash + 1);
    return path;
}


std::string replaceSuffix(const std::string& filename, const std::string& newSuffix) {
    size_t dotPos = filename.rfind(".bam");
    if (dotPos != std::string::npos) {
        return filename.substr(0, dotPos) + newSuffix;
    }
    return filename;
}

std::map<std::string, std::vector<Variant>> parseBed(const std::string& bedFile) {
    std::map<std::string, std::vector<Variant>> variantsMap;
    std::ifstream file(bedFile);
    std::string line;

    while (getline(file, line)) {
        std::istringstream iss(line);
        Variant variant;
        iss >> variant.bamFile >> variant.chr;
        std::string startPos, endPos;

        iss >> startPos >> endPos;

        if (std::isdigit(startPos[0])) {
            variant.start = std::stoll(startPos);
            variant.end = std::isdigit(endPos[0]) ? std::stoll(endPos) : variant.start;

            // Parsing either DEL/DUP or REF/ALT based on whether endPos is numeric or not
            if (std::isdigit(endPos[0])) {
                iss >> variant.varType; // For CNV
            } else {
                variant.refSeq = endPos;
                iss >> variant.altSeq; // For SNV/INDEL
                variant.varType = variant.refSeq.length() == 1 && variant.altSeq.length() == 1 ? "SNV" : "INDEL";
            }
            std::cout << " INFO: Found " + variant.varType  + " for BAM file: " << variant.bamFile << "\n";
            variantsMap[variant.bamFile].push_back(variant);
        } else {
            std::cerr << " ERROR: Malformed line in BED file: " << line << std::endl;
        }
    }
    return variantsMap;
}

std::mt19937 rng{std::random_device{}()};


int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <BAM file> <BED file> <REF FASTA> <output BAM file>" << std::endl;
        return 1;
    }

    std::string variantsFile = argv[1];
    std::string reference = argv[2];
    std::string outputDir = argv[3];

    auto variantsMap = parseBed(variantsFile);

    RefFasta ref(reference);
    std::string region;
    for (const auto& entry : variantsMap) {
        const std::string& bamFile = entry.first;
        std::cout << " INFO: Processing BAM file: " << bamFile << "\n";
        std::ifstream ifile(bamFile);
        if (!ifile) {
            std::cerr << " ERROR: BAM file does not exist: " << bamFile << std::endl;
            continue;
        }

        const std::vector<Variant>& variants = entry.second;

        BamReader reader(bamFile);
        BamRecord record;

        std::string newSuffix = ".simulated.bam";
        std::string bamName = getBasename(bamFile);
        bamName = replaceSuffix(bamName, newSuffix);
        std::string bamOut = outputDir + "/" + bamName;

        BamWriter writer(bamOut, reader.getHeader());

        for( auto& variant : variants) {
            std::cout << " INFO: Simulating variant: " << variant.varType << " at " << variant.chr << ":" << variant.start << "-" << variant.end << "\n";

            std::vector<SNV> snvs;
            if (variant.varType == "SNV" || variant.varType == "INDEL") {
                region = variant.chr + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.position - 50)) + "-" + std::to_string(variant.position + 50);
            }
            else {
                region = variant.chr + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.start)) + "-" + std::to_string(variant.end);
                snvs = identifySNVs(bamFile, ref, region, variant.chr);
            }
            reader.SetRegion(region);

            while (reader.GetNextRecord(record)) {
                if (variant.varType == "SNV") {
                    if (record.Position() <= variant.position && (record.Position() + record.Seq().length()) >= variant.position) {
                        simulateSNV(record, variant);
                        writer.WriteRecord(record);
                    }
                } 
                if (variant.varType == "INDEL") {
                    if (record.Position() <= variant.position && (record.Position() + record.Seq().length()) >= variant.position) {
                        simulateIndel(record, variant, ref);
                        writer.WriteRecord(record);
                    }
                }
                if (variant.varType == "DEL" || variant.varType == "DUP") {
                    simulateCNV(record, region, variant, snvs, ref, writer, variant.varType, 0.5);
                }
            }
        }
        BamReader fullReader(bamFile);
        BamRecord newRecord;

        while (fullReader.GetNextRecord(newRecord)) {
            bool overlapsAnyVariant = false;
            for (const auto& variant : variants) {

                std::string chr = variant.chr;
                int64_t start = variant.start;
                int64_t end = variant.end;

                if (newRecord.chrName() == chr && !(newRecord.Position() + newRecord.Seq().length() < start || newRecord.Position() > end)) {
                    overlapsAnyVariant = true;
                    break;
                }
            }
            if (!overlapsAnyVariant) {
                writer.WriteRawRecord(newRecord);
            }
        }
    }
    return 0;
}
