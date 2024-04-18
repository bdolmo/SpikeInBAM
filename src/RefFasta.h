#ifndef REFFASTA_H
#define REFFASTA_H

#include <string>
#include <htslib/faidx.h>

class RefFasta {
public:
    RefFasta(const std::string &fastaPath); // Constructor
    ~RefFasta(); // Destructor

    // Function to fetch sequence. Returns an empty string if there's an error.
    std::string fetchSequence(const std::string &chrom, int pos, int end);

private:
    faidx_t *fai;
};

#endif // REFFASTA_H

