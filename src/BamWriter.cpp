#include "BamWriter.h"
#include <stdexcept>


BamWriter::BamWriter(const std::string& filename, const bam_hdr_t* hdr) : out_(nullptr), hdr_(nullptr) {
    // Open BAM file for writing
    out_ = sam_open(filename.c_str(), "wb");
    if (!out_) {
        throw std::runtime_error("Could not open BAM file for writing: " + filename);
    }

    // Write the header
    hdr_ = bam_hdr_dup(hdr);
    if (sam_hdr_write(out_, hdr_) < 0) {
        throw std::runtime_error("Could not write header to BAM file: " + filename);
    }
    bamName = filename;
}

BamWriter::~BamWriter() {
    if (out_) sam_close(out_);
    if (hdr_) bam_hdr_destroy(hdr_);                
}

void BamWriter::WriteRecord(const BamRecord& record) {
    if (!out_) {
        throw std::runtime_error("BAM file is not open for writing");
    }

    // Convert BamRecord to bam1_t for writing
    bam1_t* b = record.ToBam1_t(); // Assuming BamRecord provides this method
    uint8_t* qual_test = bam_get_qual(b);

    if (sam_write1(out_, hdr_, b) < 0) {
        throw std::runtime_error("Could not write record to BAM file");
    }

    bam_destroy1(b); // Clean up the bam1_t struct
}

void BamWriter::WriteRawRecord(const BamRecord& record) {
    if (!out_) {
        throw std::runtime_error("BAM file is not open for writing");
    }

    // Convert BamRecord to bam1_t for writing
    bam1_t* b = record.GetBam1_t(); // Assuming BamRecord provides this method

    if (sam_write1(out_, hdr_, b) < 0) {
        throw std::runtime_error("Could not write record to BAM file");
    }

    // bam_destroy1(b); // Clean up the bam1_t struct
}
void BamWriter::CreateIndex() {
    if (sam_index_build(bamName.c_str(), 0) < 0) { // Create BAM index with default index
        throw std::runtime_error("Could not create index for BAM file: " + bamName);
    }
}