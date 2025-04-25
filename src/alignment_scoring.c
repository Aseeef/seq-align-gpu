/*
 alignment_scoring.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// Turn on debugging output by defining SEQ_ALIGN_VERBOSE
//#define SEQ_ALIGN_VERBOSE

#include <stdlib.h>
#include <stdio.h>
#include <string.h> // memset
#include <ctype.h> // tolower
#include <assert.h>
#include <stdalign.h>

#include "alignment_scoring.h"
#include "alignment_macros.h"

/**
 * Initializes a scoring_t object with the given scoring parameters.
 *
 * @param scoring          Pointer to the scoring_t structure to be initialized.
 * @param match            Score for a match between characters.
 * @param mismatch         Penalty score for a mismatch between characters.
 * @param gap_open         Penalty score for opening a gap.
 * @param gap_extend       Penalty score for extending a gap.
 * @param case_sensitive   If true, scoring is case-sensitive (uppercase and lowercase treated differently).
 */
void scoring_init(scoring_t *scoring,
                  int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool case_sensitive) {
    // Gap of length 1 has penalty (gap_open+gap_extend)
    // of length N: (gap_open + gap_extend*N)
    scoring->gap_open = gap_open;
    scoring->gap_extend = gap_extend;

    scoring->use_match_mismatch = 1;
    scoring->match = match;
    scoring->mismatch = mismatch;

    scoring->case_sensitive = case_sensitive;

    memset(scoring->swap_set, 0, sizeof(scoring->swap_set));

    scoring->min_penalty = MIN2(match, mismatch);
    scoring->max_penalty = MAX2(match, mismatch);
}

/**
 * Adds a mutation score between two specific characters.
 *
 * @param scoring          Pointer to the scoring_t structure.
 * @param a                First character involved in the mutation.
 * @param b                Second character involved in the mutation.
 * @param score            Score for the mutation (alignment between a and b).
 */
void scoring_add_mutation(scoring_t *scoring, char a, char b, int score) {
    size_t index_a = letters_to_index(a);
    size_t index_b = letters_to_index(b);
    scoring->swap_scores[index_a][index_b] = score;
    set_swap_bit(scoring, index_a, index_b);
    scoring->min_penalty = MIN2(scoring->min_penalty, score);
    scoring->max_penalty = MAX2(scoring->max_penalty, score);
}

int letters_to_index(char c) {
    if (c >= 97 && c < 123) {
        return c - 96;
    } else if (c >= 65 && c < 91) {
        return c - 64;
    } else if (c == 42) {
        return 31;
    } else {
        printf("Error: %c is not a legal character for the substitution matrix!\n", c);
        exit(1);
    }
}

char index_to_letters(int c) {
    if (c >= 1 && c < 27) {
        return c + 65;
    } else if (c == 31) {
        return '*';
    } else {
        printf("Error: %d is not a legal index for the substitution matrix!\n", c);
    }
}

/**
 * Looks up the score for aligning characters a and a batch of b's and determines if they match.
 *
 * @param scoring          Pointer to the scoring_t structure.
 * @param batch_size       The batch size
 * @param a                Query character in the alignment.
 * @param b_batch          DB batch of characters in the alignment
 * @return                 The scores for aligning a and the batch of b's.
 */
static __m256i scoring_lookup(const scoring_t *scoring, size_t batch_size, score_t a_index, score_t * b_indexes) {
    assert(batch_size == 8);

    // compute the indices we are going to use to gather
    __m256i base = _mm256_set1_epi32(a_index * 32);
    __m256i idx = _mm256_add_epi32(base, _mm256_load_si256((__m256i *) b_indexes));

    score_t * swap_scores = (int *) scoring->swap_scores;
    __m256i scores = _mm256_i32gather_epi32(swap_scores, idx, 4);

    return scores;
}

/**
 * Sets the scoring_t structure with default scoring parameters.
 *
 * @param scoring          Pointer to the scoring_t structure to be initialized with defaults.
 */
void scoring_system_default(scoring_t *scoring) {
    int match_default = 1;
    int mismatch_default = -2;
    int gap_open_default = -4;
    int gap_extend_default = -1;

    // case_sensitive = 0
    scoring_init(scoring, match_default, mismatch_default,
                 gap_open_default, gap_extend_default,
                 0);
}
