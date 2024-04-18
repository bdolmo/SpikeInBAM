#include "BamReader.h"
#include <stdexcept>

BamReader::BamReader(const std::string& bamFile): end_of_file_(false), idx_(nullptr), iter_(nullptr) {
    in_ = sam_open(bamFile.c_str(), "r");
    if (!in_) throw std::runtime_error(" ERROR: Failed to open BAM file " + bamFile);
    header_ = sam_hdr_read(in_);
    if (!header_) {
        sam_close(in_);
        throw std::runtime_error(" ERROR: Failed to read BAM header from " + bamFile);
    }
    buffer_ = bam_init1(); // Initialize the buffer

    // Load the index file for the BAM file
    idx_ = sam_index_load(in_, bamFile.c_str());
    if (!idx_) {
        throw std::runtime_error(" ERROR: Failed to load BAM index for " + bamFile);
    }
}

bam_hdr_t* BamReader::getHeader() {
    return header_;
}

BamReader::~BamReader() {
    if (header_) bam_hdr_destroy(header_);
    if (in_) sam_close(in_);
    if (buffer_) bam_destroy1(buffer_); // Make sure to free the buffer
    if (idx_) hts_idx_destroy(idx_);  // Free the index
    if (iter_) hts_itr_destroy(iter_);  // Free the iterator

}

// In BamReader.cpp
void BamReader::SetRegion(const std::string& region) {
    if (iter_) {
        hts_itr_destroy(iter_);
        iter_ = nullptr;
    }

    if (!idx_) {
        throw std::runtime_error(" ERROR: Index was not loaded properly.");
    }

    // Set the iterator for the region using sam_itr_querys
    iter_ = sam_itr_querys(idx_, header_, region.c_str());
    if (!iter_) {
        throw std::runtime_error(" ERROR: Failed to set iterator for region " + region);
    }

    // Reset EOF flag in case we're reusing this reader
    end_of_file_ = false;
}

bool BamReader::HasNext() {
    return !end_of_file_;
}

// Update GetNextRecord to work with a region
bool BamReader::GetNextRecord(BamRecord& record) {
    if (end_of_file_) {
        return false;
    }

    int ret;
    if (iter_) {  // We have a set region
        ret = sam_itr_next(in_, iter_, buffer_);
    } 
    else {  // No region set, read normally
        ret = sam_read1(in_, header_, buffer_);
    }

    if (ret >= 0) {
        record = BamRecord(buffer_, header_);
        return true;
    } 
    else {
        end_of_file_ = true;
        return false;
    }
}