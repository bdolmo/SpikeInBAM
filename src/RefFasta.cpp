#include "RefFasta.h"
#include <stdexcept>

RefFasta::RefFasta(const std::string &fastaPath) {
    // Open and index the FASTA file
    fai = fai_load(fastaPath.c_str());
    if (fai == NULL) {
        throw std::runtime_error("Could not load FASTA file: " + fastaPath);
    }
}

RefFasta::~RefFasta() {
    if (fai != NULL) {
        fai_destroy(fai); // Free the FASTA index
    }
}

std::string RefFasta::fetchSequence(const std::string &chrom, int pos, int end) {
    if (pos > end) {
        // Invalid range
        return "";
    }

    // htslib is 0-based but the end is exclusive; we adjust this accordingly
    int len; // Length of the sequence
    char *seq = faidx_fetch_seq(fai, chrom.c_str(), pos - 1, end - 1, &len);
    
    if (seq == NULL) {
        // Failed to fetch the sequence
        return "";
    }

    std::string sequence(seq);
    if (seq) free(seq); // Free memory allocated by faidx_fetch_seq

    return sequence;
}