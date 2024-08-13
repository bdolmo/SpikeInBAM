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
        std::string startPos, endPos, thirdField;

        iss >> variant.bamFile >> variant.chr >> startPos >> endPos >> thirdField;

        if (std::isdigit(startPos[0])) {
            variant.start = std::stoll(startPos);
            variant.end = variant.start; // Assume start position only, adjust based on context

            std::string copyNumberOrAltSeq;
            std::string exonName;

            if (std::isdigit(endPos[0])) {
                // This block handles CNVs
                variant.end = std::stoll(endPos);
                variant.varType = thirdField; // expecting DUP, DEL, etc.
                iss >> exonName >> copyNumberOrAltSeq;
                variant.vaf = std::stof(copyNumberOrAltSeq);

                std::cout << " INFO: Found CNV " << variant.varType << " with copy number " << variant.vaf << " for BAM file: " << variant.bamFile << "\n";
            } else {
                iss >> copyNumberOrAltSeq;

                // This block handles SNVs and INDELs
                variant.refSeq = endPos; // endPos is the reference sequence
                variant.altSeq = thirdField; // thirdField is the alternate sequence
                variant.vaf = std::stof(copyNumberOrAltSeq); // Handling VAF

                // Check VAF range
                if (variant.vaf < 0 || variant.vaf > 1) {
                    std::cerr << " ERROR: Invalid VAF value " << variant.vaf << std::endl;
                }

                // Determine variant type based on sequence lengths
                variant.varType = (variant.refSeq.length() == 1 && variant.altSeq.length() == 1) ? "SNV" : "INDEL";
                variant.end = variant.start + variant.refSeq.length() - 1; // Adjust end position for INDELs

                std::cout << " INFO: Found " << variant.varType << " with VAF " << variant.vaf << " for BAM file: " << variant.bamFile << "\n";
            }
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
    std::string suffix = argv[4];

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

        std::vector<Variant> variants = entry.second;
        std::sort(variants.begin(), variants.end(), [](const Variant& a, const Variant& b) { return a.start < b.start; });

        BamReader reader(bamFile);
        BamRecord record;

        std::string newSuffix = suffix + ".bam";
        std::string bamName = getBasename(bamFile);
        bamName = replaceSuffix(bamName, newSuffix);
        std::string bamOut = outputDir + "/" + bamName;

        BamWriter writer(bamOut, reader.getHeader());

        for (auto& variant : variants) {
            std::cout << " INFO: Simulating variant: " << variant.varType << " at " << variant.chr << ":" << variant.start << "-" << variant.end << "\n";
            std::vector<SNV> snvs;
            if (variant.varType == "SNV" || variant.varType == "INDEL") {
                region = variant.chr + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.start - 50)) + "-" + std::to_string(variant.end + 50);
            } else {
                region = variant.chr + ":" + std::to_string(std::max(static_cast<int64_t>(0), variant.start)) + "-" + std::to_string(variant.end);
                snvs = identifySNVs(bamFile, ref, region, variant.chr);
            }
            reader.SetRegion(region);

            while (reader.GetNextRecord(record)) {
                if (record.IsUnmapped()) {
                    continue;
                }
                // std::cout <<" RECORD " << record.Qname()<< std::endl;
                if (variant.varType == "SNV") {
                    if (record.Position() <= variant.start && (record.Position() + record.Seq().length()) >= variant.start) {
                        // std::cout << "here1" << std::endl;
                        // std::cout <<record.Qname() << " " <<record.Position() << " IS SNV" << std::endl;
                        simulateSNV(record, variant);
                        writer.WriteRecord(record);
                    }
                }
                if (variant.varType == "INDEL") {
                    if (record.Position() <= variant.start && (record.Position() + record.Seq().length()) >= variant.start) {
                        // std::cout << "here2" << std::endl;
                        // std::cout <<record.Qname() << " "<< record.Position() << " IS INDEL" << std::endl;
                        simulateIndel(record, variant, ref);
                        writer.WriteRecord(record);
                    }
                }
                if (variant.varType == "DEL" || variant.varType == "DUP") {
                    // std::cout << "here3" << std::endl;
                    simulateCNV(record, region, variant, snvs, ref, writer, variant.varType, 0.5);
                }
                // std::cout << "OUT" << std::endl;
            }
        }
        BamReader fullReader(bamFile);
        BamRecord newRecord;

        while (fullReader.GetNextRecord(newRecord)) {
            bool overlapsAnyVariant = false;

            auto it = std::lower_bound(variants.begin(), variants.end(), newRecord.Position(), [](const Variant& variant, int64_t pos) {
                return variant.end < pos;
            });

            while (it != variants.end() && it->start <= newRecord.Position() + newRecord.Seq().length()) {
                if (newRecord.chrName() == it->chr && !(newRecord.Position() + newRecord.Seq().length() < it->start || newRecord.Position() > it->end)) {
                    overlapsAnyVariant = true;
                    break;
                }
                ++it;
            }

            if (!overlapsAnyVariant) {
                writer.WriteRawRecord(newRecord);
            }
        }
    }
    return 0;
}
