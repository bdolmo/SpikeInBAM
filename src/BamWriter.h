#ifndef BAMWRITER_H
#define BAMWRITER_H

#include <string>
#include <htslib/sam.h>
#include "BamRecord.h" // Include your BamRecord header

class BamWriter {
public:
    BamWriter(const std::string& filename, const bam_hdr_t* hdr);
    ~BamWriter();

    void WriteRecord(const BamRecord& record);
    void WriteRawRecord(const BamRecord& record);
    bool Close(); // Add the Close function
    void CreateIndex();

private:
    std::string bamName;
    samFile* out_;
    bam_hdr_t* hdr_;
};

#endif // BAMWRITER_H
