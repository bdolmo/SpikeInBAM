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
#include "ssw_cpp.h"

struct Variant {
    std::string bamFile;
    std::string chr;
    int64_t position;

    int64_t start;      // For CNVs,start and end
    int64_t end;

    std::string refSeq;
    std::string altSeq;
    std::string varType;
    float vaf;
};

struct AlignmentResult2 {
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
    std::string compactedCigar;
    std::string extendedCigar;
};

AlignmentResult2 alignReadToContig(const std::string& read, const std::string& contig) {
    // Initialize a default scoring matrix for DNA sequences
    const int8_t match = 2;     // Match score
    const int8_t mismatch = 8;  // Mismatch penalty
    const int8_t gapOpen = 10;  // Gap open penalty
    const int8_t gapExtend = 1; // Gap extend penalty

    // Initialize the SSW aligner
    StripedSmithWaterman::Aligner aligner(match, mismatch, gapOpen, gapExtend);

    // Initialize the SSW filter (optional settings, can be adjusted as needed)
    StripedSmithWaterman::Filter filter;

    // Create an alignment object to store the result
    StripedSmithWaterman::Alignment alignment;

    // Perform the alignment
    aligner.Align(read.c_str(), contig.c_str(), contig.size(), filter, &alignment);

    // Extract the compacted CIGAR string
    std::string compactedCigar = alignment.cigar_string;

    // Generate the extended CIGAR string and adjust query_start and query_end to skip soft-clips
    std::string extendedCigar;
    int query_start = alignment.query_begin;
    int query_end = alignment.query_end;
    int ref_start = alignment.ref_begin;
    int ref_end = alignment.ref_end;

    int readPos = alignment.query_begin;
    int contigPos = alignment.ref_begin;

    // Adjust query_start by skipping leading soft-clips
    if (compactedCigar[0] == 'S') {
        size_t i = 0;
        int len = 0;
        while (isdigit(compactedCigar[i])) {
            len = len * 10 + (compactedCigar[i] - '0');
            ++i;
        }
        if (compactedCigar[i] == 'S') {
            query_start += len; // Skip the soft-clipped positions
        }
    }

    // Adjust query_end by skipping trailing soft-clips
    if (compactedCigar.back() == 'S') {
        size_t i = compactedCigar.size() - 1;
        int len = 0;
        while (isdigit(compactedCigar[i])) {
            len = len + (compactedCigar[i] - '0') * pow(10, compactedCigar.size() - i - 1);
            --i;
        }
        if (compactedCigar[i] == 'S') {
            query_end -= len; // Skip the soft-clipped positions
        }
    }

    for (size_t i = 0; i < compactedCigar.length(); ++i) {
        char op = compactedCigar[i];
        if (isdigit(op)) {
            int len = 0;
            while (isdigit(compactedCigar[i])) {
                len = len * 10 + (compactedCigar[i] - '0');
                ++i;
            }
            op = compactedCigar[i];

            for (int j = 0; j < len; ++j) {
                extendedCigar += op;
                if (op == 'M') {
                    readPos++;
                    contigPos++;
                } else if (op == 'I') {
                    readPos++;
                } else if (op == 'D') {
                    contigPos++;
                }
            }
        }
    }

    // Create the AlignmentResult struct
    AlignmentResult2 result;
    result.query_start = query_start;
    result.query_end = query_end;
    result.ref_start = ref_start;
    result.ref_end = ref_end;
    result.compactedCigar = compactedCigar;
    result.extendedCigar = extendedCigar;

    return result;
}


// std::pair<int, std::string> alignReadToContig(const std::string& read, const std::string& contig) {
//     // Initialize a default scoring matrix for DNA sequences
//     const int8_t match = 2;     // Match score
//     const int8_t mismatch = 8;  // Mismatch penalty
//     const int8_t gapOpen = 10;  // Gap open penalty
//     const int8_t gapExtend = 1; // Gap extend penalty

//     // Initialize the SSW aligner
//     StripedSmithWaterman::Aligner aligner(match, mismatch, gapOpen, gapExtend);

//     // Initialize the SSW filter (optional settings, can be adjusted as needed)
//     StripedSmithWaterman::Filter filter;

//     // Create an alignment object to store the result
//     StripedSmithWaterman::Alignment alignment;

//     // Perform the alignment
//     aligner.Align(read.c_str(), contig.c_str(), contig.size(), filter, &alignment);

//     // Convert the CIGAR vector to a string
//     std::string cigar = alignment.cigar_string;

//     // Generate the aligned sequence output using the CIGAR string
//     std::string alignedRead, alignedContig, matchLine;
//     int readPos = alignment.query_begin;
//     int contigPos = alignment.ref_begin;

//     for (size_t i = 0; i < cigar.length(); ++i) {
//         char op = cigar[i];
//         if (isdigit(op)) {
//             int len = 0;
//             while (isdigit(cigar[i])) {
//                 len = len * 10 + (cigar[i] - '0');
//                 ++i;
//             }
//             op = cigar[i];

//             if (op == 'M') { // Match/Mismatch
//                 for (int j = 0; j < len; ++j) {
//                     alignedRead += read[readPos++];
//                     alignedContig += contig[contigPos++];
//                     matchLine += (alignedRead.back() == alignedContig.back()) ? '|' : ' ';
//                 }
//             } else if (op == 'I') { // Insertion in the read
//                 for (int j = 0; j < len; ++j) {
//                     alignedRead += read[readPos++];
//                     alignedContig += '-';
//                     matchLine += ' ';
//                 }
//             } else if (op == 'D') { // Deletion in the read (insertion in the contig)
//                 for (int j = 0; j < len; ++j) {
//                     alignedRead += '-';
//                     alignedContig += contig[contigPos++];
//                     matchLine += ' ';
//                 }
//             }
//         }
//     }

//     // Calculate the total number of mismatches
//     int totalMismatches = alignment.mismatches;

//     // Print the alignment
//     std::cout << "CIGAR:   " << cigar << std::endl;
//     std::cout << "Read:    " << alignedRead << std::endl;
//     std::cout << "         " << matchLine << std::endl;
//     std::cout << "Contig:  " << alignedContig << std::endl;
//     std::cout << "Aligned length: " << alignedRead.length() << std::endl << std::endl;

//     // Return the total number of mismatches and the CIGAR string
//     return {totalMismatches, cigar};
// }



std::string str_toupper(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                // static_cast<int(*)(int)>(std::toupper)         // wrong
                // [](int c){ return std::toupper(c); }           // wrong
                // [](char c){ return std::toupper(c); }          // wrong
                   [](unsigned char c){ return std::toupper(c); } // correct
                  );
    return s;
}

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

            // writer.WriteRecord(record);
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
            // writer.WriteRecord(dupRecord);
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
    int64_t indelStart = indel.start;
    int64_t indelEnd = indel.start + indel.refSeq.length();
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

        std::random_device rd;  // Obtain a random number from hardware
        std::mt19937 gen(rd()); // Seed the generator
        std::uniform_int_distribution<> distr(0, 99); // Define the range

        int prob = indel.vaf*100;

        if (distr(gen) < prob) {
            std::string refSegment = ref.fetchSequence(record.chrName(), record.Position()-100, record.Position()+100);
            if (refSegment.empty()) {
                throw std::runtime_error("Failed to fetch reference sequence segment for alignment.");
            }
            refSegment = str_toupper(refSegment);
            // std::cout <<  modifiedSeq << " " << refSegment << std::endl;

            // AlignmentResult aln = affine_local_alignment(modifiedSeq, refSegment);
            AlignmentResult2 res = alignReadToContig(modifiedSeq, refSegment);
            // std::cout << "aln: " << aln.query_start << " " << aln.query_end << " "<< aln.ref_start << " " << aln.ref_end  << std::endl;
            // std::cout << "res: " << res.query_start << " " << res.query_end << " "<< res.ref_start << " " << res.ref_end  << std::endl;

            int num_op = 0;
            char firstOp = '\0';
            char modOp = '\0';
            std::string compact_cigar = "";

            // std::cout <<  aln.extended_cigar << std::endl;
            // std::cout <<  aln.extended_cigar << std::endl;

            for (char ntd : res.extendedCigar) {
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

            // std::cout << "COMPACT CIGAR " << compact_cigar << std::endl;

            int newRecordPosition = record.Position() - CONTEXT_SIZE + res.ref_start-1;

            record.UpdateSeq(modifiedSeq, compact_cigar);
            record.SetPosition(newRecordPosition);
        }
    }
}

void simulateSNV(BamRecord& record, const Variant& variant) {
    if (variant.refSeq.length() == 1 && variant.altSeq.length() == 1) {
        int64_t variantPosInRead = variant.start - record.Position();
        if (variantPosInRead >= 0 && variantPosInRead < static_cast<int64_t>(record.Seq().length())) {
            std::string modifiedSeq = record.Seq();
            modifiedSeq[variantPosInRead] = variant.altSeq[0];
            std::string currentCigar = record.cigarString();
            record.UpdateSeq(modifiedSeq, currentCigar);
        }
    }
}

