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
#include <x86intrin.h>

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
    int32_t *seq_a_indexes, *seq_b_batch_indexes;    // Pointers to input sequences A and B
    char *seq_a_str, **seq_b_str_batch;    // Pointers to input sequences A and B
    char *seq_a_fasta, **seq_b_fasta_batch;  // Pointers to the FASTA names
    size_t vector_size;                // the batch size of b
    size_t score_width, score_height; // Matrix dimensions: width = len(seq_a)+1, height = len(seq_b_batch[i])+1
    score_t *prev_match_scores;        // Match/mismatch array from previous row
    score_t *curr_match_scores;        // Match/mismatch array from current row
    score_t *prev_gap_a_scores;        //
    score_t *curr_gap_a_scores;        //
    score_t *prev_gap_b_scores;        //
    score_t *curr_gap_b_scores;        //
    score_t *max_scores;            // the max score of the best local alignment found
} aligner_t;

#define MATRIX_NAME(x) ((x) == MATCH ? "MATCH" : ((x) == GAP_A ? "GAP_A" : "GAP_B"))

#ifdef __cplusplus
extern "C" {
#endif

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
void aligner_update(aligner_t *aligner,
                   char *seq_a_str, char **seq_b_str_batch,
                   char * seq_a_fasta, char **seq_b_fasta_batch,
                   int32_t * seq_a_indexes, int32_t * seq_b_batch_indexes,
                   size_t len_a, size_t len_b, size_t vector_size,
                   const scoring_t *scoring);

aligner_t * aligner_create(char *seq_a_str, char **seq_b_str_batch,
                   char * seq_a_fasta, char **seq_b_fasta_batch,
                   int32_t * seq_a_indexes, int32_t * seq_b_batch_indexes,
                   size_t len_a, size_t len_b, size_t vector_size,
                   const scoring_t *scoring);

void alignment_fill_matrices(aligner_t * aligner);

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
