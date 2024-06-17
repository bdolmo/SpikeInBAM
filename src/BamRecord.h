#include <iostream>
#include <string>
#include <htslib/sam.h>
#include <vector>
#ifndef BAMRECORD_H
#define BAMRECORD_H

class BamRecord {
    public:
        BamRecord();
        BamRecord(bam1_t *b, bam_hdr_t *h);
        int chrID() const;
        int chrMateID() const;
        int64_t Position() const;
        std::string chrName() const;
        std::string chrMateName() const;
        std::string cigarString() const;
        std::string MDtag() const;
        std::string Seq() const;
        std::string Qname() const;
        int64_t matePos() const;
        // uint64_t GetID() const;

        int mapQual() const;
        void UpdateSeq(const std::string& newSeq, const std::string& newCigar);
        void SetQname(const std::string& n);
        void SetQualities(const std::string& n, int offset = 33); // Defaulting offset to Phred+33
        void SetPosition(int64_t newPos);
        // void SetCigar(const std::vector<uint32_t>& newCigar);
        std::vector<uint32_t> getCigarVector() const;

        bam1_t* ToBam1_t() const;
        bam1_t* GetBam1_t() const;
        

    private:
        uint64_t _id;
        int _chrID;
        int _chrMateID;
        int64_t _position;
        std::string _chrName;
        std::string _chrMateName;
        std::string _cigarString;
        std::string _mdString;
        std::string _seq;
        std::string _qual;
        std::string _qualString;

        int64_t _matePos;
        int _mapQual;
        int _insertSize;  // Stores the insert size (TLEN) of the record
        bam1_t *b_;
        bam_hdr_t *h_;
        int _flag;  // Stores the bitwise flags
        std::string _qname;  // Stores the query name (QNAME) of the record
};

#endif // BAMRECORD_H

