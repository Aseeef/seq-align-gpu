/*
 alignment.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Feb 2015
 */

// Turn on debugging output by defining SEQ_ALIGN_VERBOSE
//#define SEQ_ALIGN_VERBOSE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower
#include <assert.h>

#include "alignment.h"
#include "alignment_macros.h"

// Fill in traceback matrix
static void alignment_fill_matrices(aligner_t * aligner)
{
  score_t *match_scores = aligner->match_scores;
  score_t *gap_a_scores = aligner->gap_a_scores;
  score_t *gap_b_scores = aligner->gap_b_scores;
  const scoring_t *scoring = aligner->scoring;
  score_t *max_scores = aligner->max_scores;
  size_t score_width = aligner->score_width;
  size_t score_height = aligner->score_height;
  size_t i, j, b;

  int gap_open_penalty = scoring->gap_extend + scoring->gap_open;
  int gap_extend_penalty = scoring->gap_extend;

  const score_t min = 0;

  size_t seq_i, seq_j, len_i = score_width-1, len_j = score_height-1;
  size_t index, index_left, index_up, index_upleft;

  // [0][0]
  match_scores[0] = 0;
  gap_a_scores[0] = 0;
  gap_b_scores[0] = 0;

  // reset match and gap matrices
  for(i = 1; i < score_width; i++)
      match_scores[i] = gap_a_scores[i] = gap_b_scores[i] = 0;
  for(j = 1, index = score_width; j < score_height; j++, index += score_width)
      match_scores[index] = gap_a_scores[index] = gap_b_scores[index] = min;
  // reset match scores
  for (b = 0; b < aligner->b_batch_size; b++)
      max_scores[b] = min;

  // start at position [1][1]
  index_upleft = 0;
  index_up = 1;
  index_left = score_width;
  index = score_width+1;

  for (b = 0; b < aligner->b_batch_size; b++) {
      for (seq_j = 0; seq_j < len_j; seq_j++) {
          for (seq_i = 0; seq_i < len_i; seq_i++) {
              // Update match_scores[i][j] with position [i-1][j-1]
              // substitution penalty
              bool is_match;
              int substitution_penalty;

              // todo: this is a double de-reference. This was a bad idea. Optimize later.
              scoring_lookup(scoring, aligner->seq_a[seq_i], aligner->seq_b_batch[b][seq_j],
                             &substitution_penalty, &is_match);

              // substitution
              // 1) continue alignment
              // 2) close gap in seq_a
              // 3) close gap in seq_b
              score_t match_score = MAX4(match_scores[index_upleft] + substitution_penalty,
                                         gap_a_scores[index_upleft] + substitution_penalty,
                                         gap_b_scores[index_upleft] + substitution_penalty,
                                         min);
              match_scores[index] = match_score;

              // update best score
              if (max_scores[b] < match_score)
                  max_scores[b] = match_score;

              // Long arithmetic since some INTs are set to min and penalty is -ve
              // (adding as ints would cause an integer overflow)

              // Update gap_a_scores[i][j] from position [i][j-1]
              if (seq_i == len_i - 1) {
                  gap_a_scores[index]
                          = MAX4(match_scores[index_up] + gap_open_penalty,
                                 gap_a_scores[index_up] + gap_extend_penalty,
                                 gap_b_scores[index_up] + gap_open_penalty,
                                 min);
              } else
                  gap_a_scores[index] = min;

              // Update gap_b_scores[i][j] from position [i-1][j]
              if (seq_j == len_j - 1) {
                  gap_b_scores[index]
                          = MAX4(match_scores[index_left] + gap_open_penalty,
                                 gap_a_scores[index_left] + gap_open_penalty,
                                 gap_b_scores[index_left] + gap_extend_penalty,
                                 min);
              } else
                  gap_b_scores[index] = min;

              index++;
              index_left++;
              index_up++;
              index_upleft++;
          }

          index++;
          index_left++;
          index_up++;
          index_upleft++;
      }
  }
}

// Note: len_b must be same for all batches
void aligner_align(aligner_t *aligner,
                   const char *seq_a, const char **seq_b_batch,
                   size_t len_a, size_t len_b, size_t batch_size,
                   const scoring_t *scoring)
{
  aligner->scoring = scoring;
  aligner->seq_a = seq_a;
  aligner->seq_b_batch = seq_b_batch;
  aligner->b_batch_size = batch_size;
  aligner->score_width = len_a+1;
  aligner->score_height = len_b+1;

  aligner->max_scores = realloc(aligner->max_scores, sizeof(score_t) * batch_size);

  size_t new_capacity = aligner->score_width * aligner->score_height;
  if(aligner->capacity < new_capacity)
  {
    aligner->capacity = ROUNDUP2POW(new_capacity);
    size_t mem = sizeof(score_t) * aligner->capacity * batch_size;
    aligner->match_scores = realloc(aligner->match_scores, mem);
    aligner->gap_a_scores = realloc(aligner->gap_a_scores, mem);
    aligner->gap_b_scores = realloc(aligner->gap_b_scores, mem);
  }

  alignment_fill_matrices(aligner);
}

void aligner_destroy(aligner_t *aligner)
{
  if(aligner->capacity > 0) {
    free(aligner->match_scores);
    free(aligner->gap_a_scores);
    free(aligner->gap_b_scores);
  }
}


void alignment_print_matrices(const aligner_t *aligner)
{
  const score_t* match_scores = aligner->match_scores;
  const score_t* gap_a_scores = aligner->gap_a_scores;
  const score_t* gap_b_scores = aligner->gap_b_scores;

  size_t i, j;

  printf("seq_a: %.*s\nseq_b: %.*s\n",
         (int)aligner->score_width-1, aligner->seq_a,
         (int)aligner->score_height-1, aligner->seq_b);

  printf("match_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf("\t%3i", (int)ARR_LOOKUP(match_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }
  printf("gap_a_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf("\t%3i", (int)ARR_LOOKUP(gap_a_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }
  printf("gap_b_scores:\n");
  for(j = 0; j < aligner->score_height; j++)
  {
    printf("%3i:", (int)j);
    for(i = 0; i < aligner->score_width; i++)
    {
      printf("\t%3i", (int)ARR_LOOKUP(gap_b_scores, aligner->score_width, i, j));
    }
    putc('\n', stdout);
  }

  printf("match: %i mismatch: %i gapopen: %i gapexend: %i\n",
         aligner->scoring->match, aligner->scoring->mismatch,
         aligner->scoring->gap_open, aligner->scoring->gap_extend);
  printf("\n");
}


