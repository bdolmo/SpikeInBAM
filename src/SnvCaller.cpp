#include "BamReader.h"
#include "BamRecord.h"
#include "RefFasta.h"
#include <iostream>
#include <vector>
#include <map>
#include <regex>
#include <string>

struct SNV {
    std::string chr;
    int64_t position;
    std::string ref;
    std::string alt;
    int depth;
    double vaf;
    std::vector<std::string> qnameVector;
};

std::vector<SNV> identifySNVs(std::string bamFile, RefFasta& ref, const std::string& region, const std::string& chromosome) {


    // Setting the region for the BAM reader
    std::vector<SNV> snvs;

    BamReader reader(bamFile);

    BamRecord record;
    reader.SetRegion(region);

    std::map<int, std::map<char, int>> alleleCounts;
    int totalDepth = 0;

    std::map<std::string, std::vector<std::string>> qNameArray;
    while (reader.GetNextRecord(record)) {
        totalDepth++;
        std::string mdStr = record.MDtag();
        std::string seq = record.Seq();
        if (mdStr.empty()) continue;

        std::smatch matches;
        std::regex snv_regex("(\\d+)([^0-9])");
        int mdPos = 0, seqIndex = 0;

        while (std::regex_search(mdStr, matches, snv_regex)) {
            mdPos += std::stoi(matches[1]);
            seqIndex += std::stoi(matches[1]);

            if (seqIndex < seq.length()) {
                char snvBase = seq[seqIndex];
                int globalPos = record.Position() + mdPos;
                alleleCounts[globalPos][snvBase]++;
                std::string coord = chromosome + std::to_string(globalPos);
                qNameArray[coord].push_back(record.Qname());
            }
            mdPos++;
            seqIndex++;
            mdStr = matches.suffix().str(); 
        }
    }

    // Analyze allele counts to identify SNVs
    for (const auto& posMap : alleleCounts) {
        int depthAtPos = 0;       
        std::string varCoord = chromosome + ":" + std::to_string(posMap.first) + "-" + std::to_string(posMap.first);
        reader.SetRegion(varCoord);
        while (reader.GetNextRecord(record)) {
            depthAtPos++;
        }
        for (const auto& baseCount : posMap.second) {
            double vaf = static_cast<double>(baseCount.second) / depthAtPos;
            if (vaf >= 0.13 && depthAtPos >= 30) {
                std::string refBase = ref.fetchSequence(record.chrName(), posMap.first+1, posMap.first+1);
                std::string altBase(1, baseCount.first);
                SNV snv = {chromosome, posMap.first, refBase, altBase, depthAtPos, vaf};
                std::string coord = chromosome + std::to_string(posMap.first);
                for(auto& qname : qNameArray[coord]) {
                    snv.qnameVector.push_back(qname);                   
                }
                snvs.push_back(snv);
            }
        }
    }
    return snvs;
}
