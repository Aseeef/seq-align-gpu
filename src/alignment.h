/*
 alignment.h
 author: Isaac Turner <turner.isaac@gmail.com>
 url: https://github.com/noporpoise/seq-align
 May 2013
 */

#ifndef ALIGNMENT_HEADER_SEEN
#define ALIGNMENT_HEADER_SEEN

#include <string.h> // memset
#include "alignment_scoring.h"

#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) _rndup2pow64(x)
  static inline size_t _rndup2pow64(unsigned long long x) {
    // long long >=64 bits guaranteed in C99
    --x; x|=x>>1; x|=x>>2; x|=x>>4; x|=x>>8; x|=x>>16; x|=x>>32; ++x;
    return x;
  }
#endif

// Core struct for running alignments between two sequences
typedef struct
{
    const scoring_t* scoring;     // Scoring scheme used (match/mismatch/gap)
    const char *seq_a, *seq_b;    // Pointers to input sequences A and B
    size_t score_width, score_height; // Matrix dimensions: width = len(seq_a)+1, height = len(seq_b)+1
    score_t *match_scores;        // Full match/mismatch matrix (score for aligning A[i] with B[j])
    score_t *gap_a_scores;        // Matrix for gap penalties in sequence A (inserts in B)
    score_t *gap_b_scores;        // Matrix for gap penalties in sequence B (inserts in A)
    score_t max_score;            // the max score of the best local alignment found
    size_t capacity;              // Current allocated matrix size (score_width * score_height)
} aligner_t;

#define MATRIX_NAME(x) ((x) == MATCH ? "MATCH" : ((x) == GAP_A ? "GAP_A" : "GAP_B"))

#ifdef __cplusplus
extern "C" {
#endif

// Printing colour codes
extern const char align_col_mismatch[], align_col_indel[],
                  align_col_stop[];

#define aligner_init(a) (memset(a, 0, sizeof(aligner_t)))

/**
 * Runs alignment between two sequences using the given scoring parameters.
 * Populates match, gap_a, and gap_b matrices.
 *
 * Parameters:
 *   aligner   - pointer to an aligner_t struct (initialized with aligner_init)
 *   seq_a/b   - sequences to align
 *   len_a/b   - lengths of seq_a and seq_b
 *   scoring   - pointer to scoring scheme (match/mismatch/gaps)
 */
void aligner_align(aligner_t *aligner,
                   const char *seq_a, const char *seq_b,
                   size_t len_a, size_t len_b,
                   const scoring_t *scoring);

/**
 * Frees internal buffers used in the aligner.
 */
void aligner_destroy(aligner_t *aligner);

// Printing
void alignment_print_matrices(const aligner_t *aligner);

#ifdef __cplusplus
}
#endif

#endif
