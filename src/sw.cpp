#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <climits>
#include <regex>     

// Assuming match, mismatch, gap open, and gap extend scores are constants
const int MATCH = 2;
const int MISMATCH = -4;
const int GAP_OPEN = -40;
const int GAP_EXTEND = 0;

int score_match(char a, char b) {
    return a == b ? MATCH : MISMATCH;
}

int score_match_special1(char a, char b) {
    return a == b ? 0 : -24;
}

int score_match_special2(char a, char b) {
    return a == b ? -12 : -24;
}

std::string compact_cigar_string(const std::string& cigar_str) {
    if (cigar_str.empty()) return ""; // Return an empty string if the input is empty

    std::string compact_cigar;
    int count = 1; // Start counting from 1 since we always have at least one operation

    for (size_t i = 1; i <= cigar_str.size(); ++i) { // Start from the second character
        if (i < cigar_str.size() && cigar_str[i] == cigar_str[i - 1]) {
            count++; // Increment count if the current and previous characters are the same
        } else {
            // Once we reach a different character or the end of the string,
            // append the count and the operation type to the compact string
            compact_cigar += std::to_string(count) + cigar_str[i - 1];
            count = 1; // Reset count for the next operation
        }
    }
    return compact_cigar;
}


struct AlignmentResult {
    int max_i;
    int max_j;
    int max_score;
    int query_start;
    int query_end;
    int ref_start;
    int ref_end;
    std::string seq1_align;
    std::string seq2_align;
    std::string spacer;
    std::string extended_cigar;
    std::string compact_cigar;
};


AlignmentResult affine_local_alignment(const std::string& seq1, const std::string& seq2) {
    int m = seq1.length();
    int n = seq2.length();

    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> X(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> Y(m + 1, std::vector<int>(n + 1, 0));

    for (int i = 1; i <= m; ++i) {
        X[i][0] = 0;
        Y[i][0] = INT_MIN / 2; 
    }
    for (int j = 1; j <= n; ++j) {
        Y[0][j] = 0;
        X[0][j] = INT_MIN / 2;
    }

    // Filling matrices
    int max_score = 0, max_i = 0, max_j = 0;
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            M[i][j] = std::max(
                0, 
                std::max({M[i - 1][j - 1] + score_match(seq1[i - 1], seq2[j - 1]), 
                X[i - 1][j - 1]+ score_match(seq1[i - 1], seq2[j - 1]), 
                Y[i - 1][j - 1]+ score_match(seq1[i - 1], seq2[j - 1]),
                // X[i - 1][j - 1]+ score_match_special1(seq1[i - 1], seq2[j - 1]), 
                // Y[i - 1][j - 1]+ score_match_special1(seq1[i - 1], seq2[j - 1]),
            }));

            X[i][j] = std::max({
                0, 
                M[i - 1][j] + GAP_OPEN,
                X[i - 1][j] + GAP_EXTEND, 
                Y[i - 1][j]  
            });

            // Update rule for Y (deletions in seq1 or insertions in seq2)
            Y[i][j] = std::max({
                0, 
                M[i][j - 1] + GAP_OPEN,
                Y[i][j - 1] + GAP_EXTEND,
                X[i][j - 1]
            });

            // Update max score
            int current_max = std::max({M[i][j], X[i][j], Y[i][j]});
            if (current_max > max_score) {
                max_score = current_max;
                max_i = i;
                max_j = j;
            }
        }
    }
    int query_start = max_i;
    int query_end = max_i;
    int ref_start = max_j;
    int ref_end = max_j;

    // Traceback
    std::string seq1_align, seq2_align, spacer, cigar_str;
    int i = max_i, j = max_j;
    while (i > 0 || j > 0) {
        int current_max = std::max({M[i][j], X[i][j], Y[i][j]});
        if (current_max == 0) break;

        if (current_max == M[i][j]) {
            seq1_align = seq1[i - 1] + seq1_align;
            seq2_align = seq2[j - 1] + seq2_align;
            spacer = (seq1[i - 1] == seq2[j - 1] ? '|' : '*') + spacer;
            cigar_str = (seq1[i - 1] == seq2[j - 1] ? "M" : "X") + cigar_str;
            i--;
            j--;
        } else if (current_max == X[i][j]) {
            seq1_align = seq1[i - 1] + seq1_align;
            seq2_align = '-' + seq2_align;
            spacer = ' ' + spacer;
            cigar_str = "I" + cigar_str;
            i--;
        } else { 
            seq1_align = '-' + seq1_align;
            seq2_align = seq2[j - 1] + seq2_align;
            spacer = ' ' + spacer;
            cigar_str = "D" + cigar_str;
            j--;
        }
        if (i > 0) query_start = i - 1;
        if (j > 0) ref_start = j - 1;
    }

    int sum =0;
    if (i > 0) {
        sum = 1;
    }

    std::string firstSoftClip = "";
    for (int i = 0; i<query_start+sum; i++) {
        firstSoftClip += "S";
    }

    std::string secondSoftClip = "";
    for (int i = query_end;i<seq1.length();i++) {
        secondSoftClip += "S";
    }
    std::string extendedCigar = firstSoftClip + cigar_str + secondSoftClip;
    std::string compact_cigar = compact_cigar_string(cigar_str);
    // std::cout << "extended_cigar:" << extendedCigar << std::endl;
    // std::cout << seq1_align << std::endl << spacer << std::endl << seq2_align << std::endl;


    return AlignmentResult{
        max_i, max_j, max_score,
        query_start, query_end,
        ref_start, ref_end,
        seq1_align, seq2_align, spacer,
        extendedCigar,
        compact_cigar_string(extendedCigar)
    };

}

// int main() {
//     std::string seq1 = "CTACTATCGGAGAGGAAAAAGCTCGGGTATAGACCAGTATTTGTAGAGAGAACCCCCAGAGAAAAGTCGATGCGTCGAG";
//     std::string seq2 = "CTACTATCGGAGAGGAAAAAGCTCinsertionbigGAGAAAAGTCGATGCGTCGAG";
//     AlignmentResult result = affine_local_alignment(seq1, seq2);
//     std::cout << "Alignment score: " << result.max_score << std::endl;
//     std::cout << "       " << result.seq1_align << std::endl;
//     std::cout << "       " <<result.spacer << std::endl;
//     std::cout << "       " <<result.seq2_align << std::endl;
//     std::cout << "Cigar: " << result.cigar << std::endl;
//     return 0;
// }
