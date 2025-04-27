/*
 alignment_scoring.h
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

#ifndef ALIGNMENT_SCORING_HEADER_SEEN
#define ALIGNMENT_SCORING_HEADER_SEEN

#include <inttypes.h>
#include <stdbool.h>
#include <x86intrin.h>
#include <limits.h> // INT_MIN
#include <stdalign.h>

typedef int16_t score_t;
#define SCORE_MIN INT_MIN

typedef struct
{
  score_t gap_open, gap_extend;

  // If swap_score not set, should we use match/mismatch values?
  bool use_match_mismatch;
  int match, mismatch;

  bool case_sensitive;

  // Array of characters that match to everything with the same penalty (i.e. 'N's)
  uint32_t swap_set[32];
  // The penalty or the reward for a match/mismatch between two characters.
  alignas(32) int32_t swap_scores[32][32];  // swap scores int32_t for better simd

  int min_penalty, max_penalty; // min, max {match/mismatch,gapopen etc.}
} scoring_t;

#ifndef get_swap_bit
    #define get_swap_bit(scoring, a, b) \
        (((scoring)->swap_set[(size_t)(a)] >> (b)) & 1U)
    #define set_swap_bit(scoring, a, b) \
        ((scoring)->swap_set[(size_t)(a)] |= (1U << (b)))
#endif

#ifdef __cplusplus
extern "C" {
#endif

void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool case_sensitive);

void scoring_add_mutation(scoring_t* scoring, char a, char b, int score);

int letters_to_index(char c);

char index_to_letters(int c);

// Some scoring systems
void scoring_system_default(scoring_t *scoring); // DNA/RNA default

#ifdef __cplusplus
}
#endif

#endif
