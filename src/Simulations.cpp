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
#include <random>

struct Variant {
    std::string bamFile;
    std::string chr;
    int64_t position;

    int64_t start;      // For CNVs,start and end
    int64_t end;

    std::string refSeq;
    std::string altSeq;
    std::string varType;
};



void simulateCNV(BamRecord& record, const std::string& region, const Variant& variant, 
    std::vector<SNV>& snvs, RefFasta& ref, BamWriter& writer, const std::string& action, double probability) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double randomValue = dis(gen);
    bool flag = false;

    if (randomValue < 0.66) {
        flag = true;
    }

    if (randomValue < probability) {
        if (action == "DUP") {
            bool modified = false;

            std::string newQname = record.Qname() + "_dup";
            std::string modifiedSeq = record.Seq();

            writer.WriteRecord(record);
            for (auto& snv : snvs) {
                if (record.Position() <= snv.position && (record.Position() + record.Seq().length()) > snv.position) {
                    int snvPosInRead = snv.position - record.Position();
                    if (snvPosInRead < modifiedSeq.length()) {
                        modifiedSeq[snvPosInRead] = snv.alt[0];
                        modified = true;
                    }
                }
            }
        
            if (flag) {
                record.UpdateSeq(modifiedSeq, record.cigarString()); 
            }
            BamRecord dupRecord = record;
            dupRecord.SetQname(newQname);
            writer.WriteRawRecord(dupRecord);
        }
        else {
            bool modified = false;
            std::string modifiedSeq = record.Seq();
            for (auto& snv : snvs) {
                if (record.Position() <= snv.position && (record.Position() + record.Seq().length()) > snv.position) {
                    int snvPosInRead = snv.position - record.Position();
                    if (snvPosInRead < modifiedSeq.length()) {                         
                        modifiedSeq[snvPosInRead] = snv.alt[0];
                        modified = true;
                    }
                }
            }
            if (modified) {
                record.UpdateSeq(modifiedSeq, record.cigarString()); 
            }
            writer.WriteRawRecord(record);
        } 
    } 
    else {
        if (action == "DUP") {
            writer.WriteRawRecord(record);
            writer.WriteRawRecord(record);
        }
        else {
            return;
        }
    }
}

void simulateIndel(BamRecord& record, const Variant& indel, RefFasta& ref) {
    const int CONTEXT_SIZE = 100;  // Amount of context to fetch for alignment

    // using 64-bit integers since new BAM specification
    int64_t recordStart = record.Position();
    int64_t recordEnd = recordStart + record.Seq().length();
    int64_t indelStart = indel.position;
    int64_t indelEnd = indel.position + indel.refSeq.length();
    std::string modifiedSeq = record.Seq();
    bool overlapsIndel = false;

    if (record.chrName() == indel.chr && recordStart <= indelEnd && recordEnd >= indelStart) {
        int64_t overlapStart = std::max(recordStart, indelStart);
        int64_t overlapEnd = std::min(recordEnd, indelEnd);
        int64_t indelPosInRead = overlapStart - recordStart;
        int64_t indelLengthInRead = overlapEnd - overlapStart;

        if (indelStart >= recordStart && indelEnd <= recordEnd) {
            std::string newSeqPart = indel.altSeq;
            modifiedSeq.replace(indelPosInRead, indel.refSeq.length(), newSeqPart);
        }
        else {
            if (recordStart < indelStart) {
                std::string partialIndelSeq = indel.altSeq.substr(0, indelLengthInRead);
                modifiedSeq.replace(indelPosInRead, indelLengthInRead, partialIndelSeq);
            }
            if (recordEnd > indelEnd) {
                std::string partialIndelSeq = indel.altSeq.substr(0, indelLengthInRead); 
                modifiedSeq.replace(indelPosInRead, indelLengthInRead, partialIndelSeq);
            }
        }
        overlapsIndel = true;
    }
    if (overlapsIndel) {
        std::string refSegment = ref.fetchSequence(record.chrName(), record.Position()-100, record.Position()+100);
        if (refSegment.empty()) {
            throw std::runtime_error("Failed to fetch reference sequence segment for alignment.");
        }
        AlignmentResult aln = affine_local_alignment(modifiedSeq, refSegment);

        int num_op = 0;
        char firstOp = '\0';
        char modOp = '\0';
        std::string compact_cigar = "";

        for (char ntd : aln.extended_cigar) {
            if (firstOp == '\0') {
                firstOp = ntd;
                modOp = ntd;
                num_op = 1;
            } 
            else if (firstOp == ntd || (firstOp == 'M' && ntd == 'X') || (firstOp == 'X' && ntd == 'M')) {
                num_op++;
                if ((firstOp == 'M' && ntd == 'X') || (firstOp == 'X' && ntd == 'M')) {
                    modOp = 'M';
                }
            } 
            else {
                compact_cigar += std::to_string(num_op) + modOp;
                firstOp = ntd;
                modOp = ntd;
                num_op = 1;
            }
        }
        if (num_op > 0) {
            compact_cigar += std::to_string(num_op) + modOp;
        }
        int newRecordPosition = record.Position() - CONTEXT_SIZE + aln.ref_start;
        record.UpdateSeq(modifiedSeq, compact_cigar);
        record.SetPosition(newRecordPosition);
    }
}

void simulateSNV(BamRecord& record, const Variant& variant) {
    if (variant.refSeq.length() == 1 && variant.altSeq.length() == 1) {
        int64_t variantPosInRead = variant.position - record.Position();
        if (variantPosInRead >= 0 && variantPosInRead < static_cast<int64_t>(record.Seq().length())) {
            std::string modifiedSeq = record.Seq();
            modifiedSeq[variantPosInRead] = variant.altSeq[0];
            std::string currentCigar = record.cigarString();
            record.UpdateSeq(modifiedSeq, currentCigar);
        }
    }
}

