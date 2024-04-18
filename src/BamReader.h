#ifndef BAMREADER_H
#define BAMREADER_H

#include <string>
#include <htslib/sam.h>
#include "BamRecord.h"

class BamReader {
public:
    explicit BamReader(const std::string& filename);
    ~BamReader();
    bool GetNextRecord(BamRecord& record); 
    bool HasNext(); 
    void SetRegion(const std::string& region);  // Add this line
    bam_hdr_t* getHeader();
private:
    samFile *in_;
    bam_hdr_t *header_;
    bam1_t *buffer_;
    bool end_of_file_;
    hts_idx_t *idx_;
    hts_itr_t *iter_;
    
};

#endif
