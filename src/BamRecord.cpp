#include "BamRecord.h"
#include <stdexcept>
#include <string>
#include <cstring>  // For memcpy and memset
#include <cstdlib>  // For malloc, free
#include <stdexcept> // For std::runtime_error
#include <htslib/sam.h>
#include <regex>
BamRecord::BamRecord() {
}

BamRecord::BamRecord(bam1_t *b, bam_hdr_t *h) : b_(b), h_(h) {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer received for BAM record");
    }
    if (h_ == nullptr) {
        throw std::runtime_error("Null pointer received from BAM header");
    }

    _position = b->core.pos;
    _chrID = b->core.tid;
    _chrMateID = b->core.mtid;

    if (b->core.tid >= 0) {
        _chrName = std::string(h->target_name[b->core.tid]);
    } 
    else {
        _chrName = "*";
    }
    if (b->core.mtid >= 0) {
        _chrMateName = std::string(h->target_name[b->core.mtid]);
    } 
    else {
        _chrMateName = "*";
    }
    _qname = bam_get_qname(b);  // Extract QNAME from the bam1_t structure
    _flag = b->core.flag; // Derive the flag from the bam1_t structure
    const uint32_t *cigar = bam_get_cigar(b);
    for (unsigned i = 0; i < b->core.n_cigar; ++i) {
        int cigarLen = bam_cigar_oplen(cigar[i]);
        char cigarOp = bam_cigar_opchr(cigar[i]);
        _cigarString += std::to_string(cigarLen) + cigarOp;
    }
    _insertSize = b->core.isize;
    // Borrowed from SeqLib
    uint8_t* p = bam_aux_get(b, "MD");
    if (!p) {
        _mdString = "";
    }
    else {
        char* pp = bam_aux2Z(p);
        if (!pp) {
            _mdString = "";
        }
        else {
            _mdString = std::string(pp);
        }       
    }

    _mapQual = b->core.qual;
    _matePos = b->core.mpos;

    uint8_t *seq = bam_get_seq(b); // Get compressed sequence
    int32_t len = b->core.l_qseq; // Get length of the sequence
    _seq.reserve(len); // Reserve space to avoid multiple reallocations

    for (int i = 0; i < len; ++i) {
        uint8_t base = bam_seqi(seq, i); // Get base in 4-bit integer format
        switch (base) {
            case 1: _seq += 'A'; break;
            case 2: _seq += 'C'; break;
            case 4: _seq += 'G'; break;
            case 8: _seq += 'T'; break;
            case 15: _seq += 'N'; break; // 15 represents 'N' in BAM format
            default: _seq += '?'; break; // Just in case there is an unknown base
        }
    }

    // Set the quality
    uint8_t* _qual = bam_get_qual(b);

    // std::string _qualString;
    _qualString.reserve(len);
    for (int i = 0; i < len; ++i) {
        // Convert the Phred quality score to the ASCII representation (Phred+33)
        _qualString.push_back(_qual[i]);
    }
}

std::vector<uint32_t> BamRecord::getCigarVector() const {
    std::vector<uint32_t> cigarVector;
    std::regex cigarRegex("(\\d+)([MIDNSHP=X])"); // Regular expression to parse CIGAR string
    std::smatch match;
    std::string tempCigar = _cigarString;

    // Use regular expression to extract each CIGAR operation
    while (std::regex_search(tempCigar, match, cigarRegex)) {
        uint32_t len = std::stoi(match[1].str()); // Length of the CIGAR operation
        char opChar = match[2].str()[0]; // Character representing the CIGAR operation
        uint32_t op; // Numeric code for the CIGAR operation

        // Convert the CIGAR operation character to BAM CIGAR operation code
        switch (opChar) {
            case 'M': op = BAM_CMATCH; break;
            case 'I': op = BAM_CINS; break;
            case 'D': op = BAM_CDEL; break;
            case 'N': op = BAM_CREF_SKIP; break;
            case 'S': op = BAM_CSOFT_CLIP; break;
            case 'H': op = BAM_CHARD_CLIP; break;
            case 'P': op = BAM_CPAD; break;
            case '=': op = BAM_CEQUAL; break;
            case 'X': op = BAM_CDIFF; break;
            default: throw std::runtime_error("Unknown CIGAR operation: " + std::string(1, opChar));
        }

        // Encode the CIGAR operation into the 32-bit format used by BAM files
        uint32_t encodedCigar = (len << BAM_CIGAR_SHIFT) | op;
        cigarVector.push_back(encodedCigar);

        // Continue searching from the end of the last match
        tempCigar = match.suffix().str();
    }

    return cigarVector;
}
// Helper function for converting base characters to their 4-bit encoded values
uint8_t base2val(char base) {
    switch (base) {
        case 'A': case 'a': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 4;
        case 'T': case 't': return 8;
        default: return 15; // 'N' or unknown base
    }
}

bam1_t* BamRecord::GetBam1_t() const {
    return b_;
}

bam1_t* BamRecord::ToBam1_t() const {
    // Allocate a new bam1_t structure

    // return b_;

    bam1_t* b_new = bam_init1();

    if (!b_new) throw std::runtime_error("Failed to allocate bam1_t structure");

    // Convert the CIGAR string to a vector of uint32_t
    std::vector<uint32_t> cigarVector = getCigarVector();
    size_t qname_len = _qname.length() + 1; // +1 for null terminator
    size_t seq_len = (_seq.length() + 1) / 2; // Sequence is 4-bit encoded
    size_t qual_len = _seq.length(); // Quality score length equals sequence length
    size_t cigar_bytes = cigarVector.size() * 4; // 4 bytes per CIGAR op

    // Calculate total data size needed
    size_t total_data_size = qname_len + cigar_bytes + seq_len + qual_len;

    // Ensure space for all data components is allocated
    if (b_new->m_data < total_data_size) {
        b_new->data = (uint8_t*)realloc(b_new->data, total_data_size);
        if (!b_new->data) {
            bam_destroy1(b_new); // Free memory if reallocation failed
            throw std::runtime_error("Failed to reallocate memory for bam1_t structure");
        }
        b_new->m_data = total_data_size;
    }

    // Set basic fields from BamRecord to bam1_t
    b_new->core.tid = _chrID;
    b_new->core.pos = _position;
    b_new->core.qual = _mapQual;
    b_new->core.mtid = _chrMateID;
    b_new->core.l_qseq = _seq.length();
    b_new->core.n_cigar = cigarVector.size();
    b_new->core.flag = _flag;
    b_new->l_data = total_data_size;
    b_new->core.isize = _insertSize;
    b_new->core.mpos = _matePos;
    b_new->core.l_qname = qname_len;

    if (!_mdString.empty()) {
        // Append the MD tag as a null-terminated string
        if (bam_aux_append(b_new, "MD", 'Z', _mdString.length() + 1, 
                        reinterpret_cast<const uint8_t*>(_mdString.c_str())) < 0) {
            throw std::runtime_error("Failed to append MD tag to bam1_t structure");
        }
    }

    // Set the query name (QNAME)
    if (qname_len > 0) {
        std::memcpy(b_new->data, _qname.c_str(), qname_len);
    }

    // Set the CIGAR operations
    uint32_t* cigar_ptr = bam_get_cigar(b_new);
    for (size_t i = 0; i < cigarVector.size(); ++i) {
        cigar_ptr[i] = cigarVector[i];
    }

    // Set the sequence
    uint8_t* seq = bam_get_seq(b_new);
    for (size_t i = 0; i < _seq.length(); ++i) {
        seq[i >> 1] = (i & 1) ? (seq[i >> 1] | base2val(_seq[i])) : (base2val(_seq[i]) << 4);
    }

    // // Set the quality
    uint8_t* quality = bam_get_qual(b_new);
    for (size_t i = 0; i < qual_len; ++i) {
        // quality[i] = quality_test[i];
        quality[i] =  _qualString[i];
    }

    // free(q);
    memcpy(bam_get_qual(b_new), quality, qual_len); // Copy converted quality scores back into the BAM structure


    return b_new;
}

void BamRecord::SetPosition(int64_t newPos) {
    if (b_ == nullptr) {
        throw std::runtime_error("Null pointer encountered in BamRecord");
    }
    _position = newPos; // Update the internal representation
    b_->core.pos = newPos; // Update the position in the bam1_t struct
}


void BamRecord::UpdateSeq(const std::string& seq, const std::string& cigar) {

    // Borrowed from SeqLib
    int new_size = b_->l_data - ((b_->core.l_qseq+1)>>1) - b_->core.l_qseq + ((seq.length()+1)>>1) + seq.length();    
    int old_aux_spot = (b_->core.n_cigar<<2) + b_->core.l_qname + ((b_->core.l_qseq + 1)>>1) + b_->core.l_qseq;
    int old_aux_len = bam_get_l_aux(b_); //(b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;

    // copy out all the old data
    uint8_t* oldd = (uint8_t*)malloc(b_->l_data);
    memcpy(oldd, b_->data, b_->l_data);

    // clear out the old data and alloc the new amount
    free(b_->data);
    b_->data = (uint8_t*)calloc(new_size, sizeof(uint8_t)); 

    // add back the qname and cigar
    memcpy(b_->data, oldd, b_->core.l_qname + (b_->core.n_cigar<<2));

    // update the sizes
    // >>1 shift is because only 4 bits needed per ATCGN base
    b_->l_data = new_size;
    b_->core.l_qseq = seq.length();

    // allocate the sequence
    uint8_t* m_bases = b_->data + b_->core.l_qname + (b_->core.n_cigar<<2);
    int slen = seq.length();

    _seq = "";
    _seq.reserve(slen);
    for (int i = 0; i < slen; ++i) {

        uint8_t base = 15;
        if (seq.at(i) == 'A') {
        base = 1;
        _seq+= 'A';
        }
        else if (seq.at(i) == 'C') {
        base = 2;
        _seq+= 'C';
        }
        else if (seq.at(i) == 'G') {
        base = 4;
        _seq+= 'G';
        }
        else if (seq.at(i) == 'T') {
        base = 8;
        _seq+= 'T';
        }
        else {
        _seq+= 'N';
        }
        
        m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
        m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding
    }
    _cigarString = cigar;
    std::vector<uint32_t> cigarVector = getCigarVector();
    size_t cigar_bytes = cigarVector.size() * 4; // 4 bytes per CIGAR op

    // Set the CIGAR operations
    uint32_t* cigar_ptr = bam_get_cigar(b_);
    for (size_t i = 0; i < cigarVector.size(); ++i) {
        cigar_ptr[i] = cigarVector[i];
    }

    uint8_t* quality = bam_get_qual(b_);

    for (size_t i = 0; i < slen; ++i) {
        quality[i] = _qualString[i];
    }

    // free(q);
    memcpy(bam_get_qual(b_), quality, slen); // Copy converted quality scores back into the BAM structure


    // add the aux data
    uint8_t* t = bam_get_aux(b_);
    memcpy(t, oldd + old_aux_spot, old_aux_len);

    // reset the max size
    b_->m_data = b_->l_data;

    free(oldd); //just added
}

void BamRecord::SetQname(const std::string& n) {
    // Copy out the non-qname data
    size_t nonq_len = b_->l_data - b_->core.l_qname;
    uint8_t* nonq = (uint8_t*)malloc(nonq_len);
    if (!nonq) throw std::runtime_error("Failed to allocate memory for non-qname data");
    memcpy(nonq, b_->data + b_->core.l_qname, nonq_len);

    // Clear the old data and allocate the new amount
    free(b_->data);
    b_->data = (uint8_t*)calloc(nonq_len + n.length() + 1, 1);
    if (!b_->data) {
        free(nonq); // Avoid memory leak
        throw std::runtime_error("Failed to allocate memory for new data");
    }
    
    // Add in the new qname
    memcpy(b_->data, (uint8_t*)n.c_str(), n.length() + 1); // +1 for '\0' terminator

    // Update the sizes
    b_->l_data = nonq_len + n.length() + 1;
    b_->core.l_qname = n.length() + 1;
    
    // Copy over the old data after the new qname
    memcpy(b_->data + b_->core.l_qname, nonq, nonq_len);
    free(nonq); // Free the temporary buffer

    // Reset the max size
    b_->m_data = b_->l_data;
}

void BamRecord::SetQualities(const std::string& n, int offset) {
    if (!n.empty() && n.length() != static_cast<size_t>(b_->core.l_qseq))
        throw std::invalid_argument("New quality score string should be the same length as sequence length");
    
    // Length of qual is always same as seq. If empty qual, just set first bit of qual to 0 (equivalent to '!')
    if (n.empty()) {
        uint8_t* r = bam_get_qual(b_);
        r[0] = 0; // In BAM format, a quality score of 0 is used to represent '!'
        return;
    }

    // Convert quality string to numeric values
    char * q = strdup(n.c_str());
    std::cout << n << std::endl;
    if (!q) throw std::runtime_error("Failed to duplicate quality string");
    for (size_t i = 0; i < n.length(); ++i) {
        q[i]; // Adjust quality scores based on the offset (typically 33 or 64)
    }
    memcpy(bam_get_qual(b_), q, n.length()); // Copy converted quality scores back into the BAM structure
    free(q); // Free the duplicated quality string
}

int BamRecord::chrID() const {
    return _chrID;
}

int BamRecord::chrMateID() const {
    return _chrMateID;
}

int64_t BamRecord::Position() const{
    return _position;
}

std::string BamRecord::chrName() const {
    return _chrName;
}

std::string BamRecord::chrMateName() const {
    return _chrMateName;
}

std::string BamRecord::cigarString() const {
    return _cigarString;
}

int BamRecord::mapQual() const {
    return _mapQual;
}

std::string BamRecord::MDtag() const {
    return _mdString;
}

std::string BamRecord::Seq() const {
    return _seq;
}

std::string BamRecord::Qname() const {
    // Check if the bam1_t structure is valid
    if (!b_ || !(b_->data)) {
        throw std::runtime_error("Invalid BAM record");
    }

    // The QNAME is stored at the beginning of the data section of bam1_t structure
    // and is NULL-terminated.
    const char* qname = bam_get_qname(b_);

    // Convert C string to C++ string and return
    return std::string(qname);
}


// uint64_t BamRecord::GetID() const {
//     return _id;
// }
