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

typedef int32_t score_t;
#define SCORE_MIN INT_MIN

typedef struct
{
  int gap_open, gap_extend;

  // If swap_score not set, should we use match/mismatch values?
  bool use_match_mismatch;
  int match, mismatch;

  bool case_sensitive;

  // Array of characters that match to everything with the same penalty (i.e. 'N's)
  uint32_t swap_set[256][256/32];
  // The penalty or the reward for a match/mismatch between two characters.
  score_t swap_scores[256][256];

  int min_penalty, max_penalty; // min, max {match/mismatch,gapopen etc.}
} scoring_t;

#ifndef bitset32_get
  #define bitset32_get(arr,idx)   (((arr)[(idx)>>5] >> ((idx)&31)) & 0x1)
  #define bitset32_set(arr,idx)   ((arr)[(idx)>>5] |=   (1<<((idx)&31)))
  #define bitset32_clear(arr,idx) ((arr)[(idx)>>5] &=  ~(1<<((idx)&31)))
#endif

#define get_swap_bit(scoring,a,b) bitset32_get((scoring)->swap_set[(size_t)(a)],b)
#define set_swap_bit(scoring,a,b) bitset32_set((scoring)->swap_set[(size_t)(a)],b)

#ifdef __cplusplus
extern "C" {
#endif

void scoring_init(scoring_t* scoring, int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool case_sensitive);

void scoring_add_mutation(scoring_t* scoring, char a, char b, int score);

__m256i scoring_lookup(const scoring_t* scoring, size_t batch_size, char a, char * b);

// Some scoring systems
void scoring_system_default(scoring_t *scoring); // DNA/RNA default

#ifdef __cplusplus
}
#endif

#endif
