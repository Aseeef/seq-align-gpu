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
#include <assert.h>
#include <omp.h>
#include <x86intrin.h> // <immintrin.h> doesnt have gather

#include "alignment.h"
#include "alignment_macros.h"

/**
 * Looks up the score for aligning characters a and a batch of b's and determines if they match.
 *
 * @param scoring          Pointer to the scoring_t structure.
 * @param batch_size       The batch size
 * @param a                Query character in the alignment.
 * @param b_batch          DB batch of characters in the alignment
 * @return                 The scores for aligning a and the batch of b's.
 */
inline static __m256i scoring_lookup(const scoring_t *scoring, int32_t a_index, int32_t *b_indexes) {
    // Load the base address (the row) of the swap_scores for the given `a_index`
    const int32_t *swap_scores = scoring->swap_scores[a_index];

    // Load the indices for the batch
    // must load two since simd batch is 16, and
    // simd gather cant be used with shorts
    __m256i idx_low = _mm256_load_si256((__m256i *)b_indexes);
    __m256i idx_high = _mm256_load_si256((__m256i *)(b_indexes + 8));

    // gather the scores from the swap_scores array
    __m256i scores_low = _mm256_i32gather_epi32(swap_scores, idx_low, 4);
    __m256i scores_high = _mm256_i32gather_epi32(swap_scores, idx_high, 4);

    // clamps the scores into a vectors of shorts
    __m256i packed_scores_16 = _mm256_packs_epi32(scores_low, scores_high);

    return packed_scores_16;
}

// Fill in traceback matrix for an ENTIRE BATCH
static void alignment_fill_matrices(aligner_t * aligner)
{
  score_t *match_scores = aligner->match_scores;
  score_t *gap_a_scores = aligner->gap_a_scores;
  score_t *gap_b_scores = aligner->gap_b_scores;
  const scoring_t *scoring = aligner->scoring;
  size_t score_width = aligner->score_width;
  size_t score_height = aligner->score_height;
  size_t i, j;

  __m256i gap_open_penalty = _mm256_set1_epi16(scoring->gap_extend + scoring->gap_open);
  __m256i gap_extend_penalty = _mm256_set1_epi16(scoring->gap_extend);
  __m256i min_v = _mm256_setzero_si256();

  const size_t batch_size = aligner->b_batch_size;

  // null checks
  assert(match_scores != NULL);
  assert(gap_a_scores != NULL);
  assert(gap_b_scores != NULL);
  assert(scoring != NULL);

  size_t seq_i, seq_j, len_i = score_width-1, len_j = score_height-1;
  size_t index, index_left, index_right, index_up, index_upleft;

  // reset match scores
  __m256i max_scores_vec = _mm256_setzero_si256();

  // reset match and gap matrices
  // The shape is height x width x b
  // I need to reset first col and first row of each batch
  __m256i zero_v = _mm256_setzero_si256();
  for (i = 0; i < score_width; i++) {
      size_t offset = i * batch_size;
      _mm256_store_si256((__m256i *)(match_scores + offset), zero_v);
      _mm256_store_si256((__m256i *)(gap_a_scores + offset), zero_v);
      _mm256_store_si256((__m256i *)(gap_b_scores + offset), zero_v);
  }
  for (j = 0; j < score_height; j++) {
      size_t offset = j * score_width * batch_size;
      _mm256_store_si256((__m256i *)(match_scores + offset), zero_v);
      _mm256_store_si256((__m256i *)(gap_a_scores + offset), zero_v);
      _mm256_store_si256((__m256i *)(gap_b_scores + offset), zero_v);
  }

  // start at position [1][1]
  index_upleft = 0;
  index_up = batch_size;
  index_left = score_width*batch_size;
  index = (score_width*batch_size)+batch_size;
  index_right = (score_width*batch_size)+(2*batch_size);

  // todo: Refactor the loops and combine with OMP to do parallel anti-diagonal wavefront
  //  calculation. First step with the OMP would be to refactor using naive wavefront. Then combine it with
  //  wavefront + possibly blocking using strips (limiting the number of rows being operated on for cache friendliness
  //  since).
  //  Actually one thing I have to remember is vectors are 256-bits (32 bytes), and cache line is 512 bits (64 bytes).
  //  So anytime I hit cache I'm actually grabbing two vectors.
  //  Number of rows should be AT LEAST as must as number of cores to fully utilize CPU. Works out for me since
  //  CPUs with larger core counts tend to have larger cache. This is not as good as Reza's implementation since he split
  //  the cores across different database entries resulting in better L1/L2 cache utilization. Another thing I need to
  //  watch out for with OpenMP is false sharing. Something to think about later.

  size_t next_a_index;
  // int index = (h * width + w) * batch_size;
  for (seq_j = 0; seq_j < len_j; seq_j++) {

      for (seq_i = 0; seq_i < len_i; seq_i++) {

          if (seq_i + 1 < len_i) {
              next_a_index = aligner->seq_a_indexes[seq_i + 1];
          } else if (seq_j + 1 < len_j) {
              next_a_index = aligner->seq_a_indexes[0];
          }
          // Calculate the address of the start of the 128-byte block for the future a_index
          // this prefetches into L1 cache
          char const* prefetch_addr = (char const*)(scoring->swap_scores + next_a_index * 32);
          // I need the next 128 bytes for the next loop of the scoring lookup
          _mm_prefetch(prefetch_addr, _MM_HINT_T1);
          _mm_prefetch(prefetch_addr + 64, _MM_HINT_T1);

          // substitution penalty
          __m256i substitution_penalty = scoring_lookup(scoring, aligner->seq_a_indexes[seq_i], aligner->seq_b_batch_indexes + (seq_j*batch_size));

          // Prefetch gap_a_scores, gap_b_scores, match_scores at the upcoming index
          char const* gap_a_scores_ptr = (char const*) (gap_a_scores + index_right);
          char const* gap_b_scores_ptr = (char const*) (gap_b_scores + index_right);
          char const* match_scores_ptr = (char const*) (match_scores + index_right);
          _mm_prefetch(gap_a_scores_ptr, _MM_HINT_T0);
          _mm_prefetch(gap_b_scores_ptr, _MM_HINT_T0);
          _mm_prefetch(match_scores_ptr, _MM_HINT_T0);

          // 1) continue alignment
          // 2) close gap in seq_a
          // 3) close gap in seq_b

          // Update match_scores[i][j]
//          score_t match_score = MAX4(match_scores[index_upleft] + substitution_penalty,
//                                     gap_a_scores[index_upleft] + substitution_penalty,
//                                     gap_b_scores[index_upleft] + substitution_penalty,
//                                     min);

          __m256i match_score = _mm256_load_si256((__m256i *)(match_scores + index_upleft));
          match_score = _mm256_add_epi16(match_score, substitution_penalty);

          __m256i gap_a_score = _mm256_load_si256((__m256i *)(gap_a_scores + index_upleft));
          gap_a_score = _mm256_add_epi16(gap_a_score, substitution_penalty);
          __m256i gap_b_score = _mm256_load_si256((__m256i *)(gap_b_scores + index_upleft));
          gap_b_score = _mm256_add_epi16(gap_b_score, substitution_penalty);

          __m256i max_ab = _mm256_max_epi16(gap_a_score, gap_b_score);
          __m256i max_mab = _mm256_max_epi16(match_score, max_ab);
          match_score = _mm256_max_epi16(max_mab, min_v);

          // save
          _mm256_store_si256((__m256i *)(match_scores + index), match_score);

          // update best score
          // equal to: max_scores_vec[i] = (match_score[i] > max_scores_vec[i]) ? match_score[i] : max_scores_vec[i];
          __m256i mask = _mm256_cmpgt_epi16(match_score, max_scores_vec);
          max_scores_vec = _mm256_blendv_epi8(max_scores_vec, match_score, mask);

          // Long arithmetic since some INTs are set to min and penalty is -ve
          // (adding as ints would cause an integer overflow)

          // Update gap_a_scores[i][j]
//          gap_a_scores[index]
//                  = MAX4(match_scores[index_up] + gap_open_penalty,
//                         gap_a_scores[index_up] + gap_extend_penalty,
//                         gap_b_scores[index_up] + gap_open_penalty,
//                         min);
          match_score = _mm256_load_si256((__m256i *)(match_scores + index_up));
          match_score = _mm256_add_epi16(match_score, gap_open_penalty);
          gap_a_score = _mm256_load_si256((__m256i *)(gap_a_scores + index_up));
          gap_a_score = _mm256_add_epi16(gap_a_score, gap_extend_penalty);
          gap_b_score = _mm256_load_si256((__m256i *)(gap_b_scores + index_up));
          gap_b_score = _mm256_add_epi16(gap_b_score, gap_open_penalty);
          gap_a_score = _mm256_max_epi16(gap_a_score, match_score);
          gap_a_score = _mm256_max_epi16(gap_a_score, gap_b_score);
          gap_a_score = _mm256_max_epi16(gap_a_score, min_v);
          _mm256_store_si256((__m256i *)(gap_a_scores + index), gap_a_score);

          // Update gap_b_scores[i][j]
//          gap_b_scores[index]
//                  = MAX4(match_scores[index_left] + gap_open_penalty,
//                         gap_a_scores[index_left] + gap_open_penalty,
//                         gap_b_scores[index_left] + gap_extend_penalty,
//                         min);
          match_score = _mm256_load_si256((__m256i *)(match_scores + index_left));
          match_score = _mm256_add_epi16(match_score, gap_open_penalty);
          gap_a_score = _mm256_load_si256((__m256i *)(gap_a_scores + index_left));
          gap_a_score = _mm256_add_epi16(gap_a_score, gap_open_penalty);
          gap_b_score = _mm256_load_si256((__m256i *)(gap_b_scores + index_left));
          gap_b_score = _mm256_add_epi16(gap_b_score, gap_extend_penalty);
          gap_b_score = _mm256_max_epi16(gap_b_score, match_score);
          gap_b_score = _mm256_max_epi16(gap_b_score, gap_a_score);
          gap_b_score = _mm256_max_epi16(gap_b_score, min_v);
          _mm256_store_si256((__m256i *)(gap_b_scores + index), gap_b_score);

          index += batch_size;
          index_left += batch_size;
          index_up += batch_size;
          index_upleft += batch_size;
          index_right += batch_size;
      }
      // need to increment again to make everyone
      // get to the next row
      index += batch_size;
      index_left += batch_size;
      index_up += batch_size;
      index_upleft += batch_size;
      index_right += batch_size;
  }

    //printf("Matrix Index [3][2] = %d\n", match_scores[ARR_3D_index(score_width, batch_size, 3, 2, 7)]);

  // put back the max scores in this batch
  assert(aligner->match_scores != NULL);
  _mm256_storeu_si256((__m256i *)(aligner->max_scores), max_scores_vec);
}

// Note: len_b must be same for all batches
void aligner_align(aligner_t *aligner,
                   char *seq_a, char **seq_b_batch,
                   score_t * seq_a_indexes, score_t * seq_b_batch_indexes,
                   size_t len_a, size_t len_b, size_t batch_size,
                   const scoring_t *scoring)
{
  aligner->scoring = scoring;
  aligner->seq_a_str = seq_a;
  aligner->seq_b_str_batch = seq_b_batch;
  aligner->seq_a_indexes = seq_a_indexes;
  aligner->seq_b_batch_indexes = seq_b_batch_indexes;
  aligner->b_batch_size = batch_size;
  aligner->score_width = len_a+1; // for col of all zeros
  aligner->score_height = len_b+1; // for the row of all zeros

  // The first allocation is expected to be the largest width*height.
  assert(aligner->capacity == 0 || ROUNDUP2POW(aligner->score_width * aligner->score_height) <= aligner->capacity);

  if(aligner->max_scores == NULL)
  {
    size_t capacity = ROUNDUP2POW(aligner->score_width * aligner->score_height);
    aligner->capacity = capacity;
    aligner->max_scores = aligned_alloc(32, sizeof(score_t) * batch_size);
    size_t mem = sizeof(score_t) * capacity * batch_size;
    aligner->match_scores = aligned_alloc(32, mem);
    aligner->gap_a_scores = aligned_alloc(32, mem);
    aligner->gap_b_scores = aligned_alloc(32, mem);
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


void alignment_print_matrices(const aligner_t *aligner, size_t batch_size)
{
  const score_t* match_scores = aligner->match_scores;
  const score_t* gap_a_scores = aligner->gap_a_scores;
  const score_t* gap_b_scores = aligner->gap_b_scores;

  size_t i, j, b;


    // Note: Yes, the indexing here is not cache friendly but it's fine, this code is only
  // used to debug. So I rather not refactor it...
  for (b = 0 ; b < batch_size; b++) {

      printf("(batch no: %zi/%zi, seq_a: %.*s\nseq_b: %.*s\n",
             b + 1, batch_size,
             (int)aligner->score_width-1, aligner->seq_a_str,
             strlen(aligner->seq_b_str_batch[b]), aligner->seq_b_str_batch[b]);

      printf("\n");

      printf("match_scores:\n");
      for(j = 0; j < aligner->score_height; j++)
      {
          printf("%3i:", (int)j);
          for(i = 0; i < aligner->score_width; i++)
          {
              printf("\t%3i", (int) match_scores[ARR_3D_index(aligner->score_width, batch_size, j, i, b)]);
              //printf("\t%3i", (int)ARR_LOOKUP(match_scores, aligner->score_width, i, j));
          }
          putc('\n', stdout);
      }
      printf("gap_a_scores:\n");
      for(j = 0; j < aligner->score_height; j++)
      {
          printf("%3i:", (int)j);
          for(i = 0; i < aligner->score_width; i++)
          {
              printf("\t%3i", (int) gap_a_scores[ARR_3D_index(aligner->score_width, batch_size, j, i, b)]);
          }
          putc('\n', stdout);
      }
      printf("gap_b_scores:\n");
      for(j = 0; j < aligner->score_height; j++)
      {
          printf("%3i:", (int)j);
          for(i = 0; i < aligner->score_width; i++)
          {
              printf("\t%3i", (int) gap_b_scores[ARR_3D_index(aligner->score_width, batch_size, j, i, b)]);
          }
          putc('\n', stdout);
      }

      printf("match: %i mismatch: %i gapopen: %i gapexend: %i\n",
             aligner->scoring->match, aligner->scoring->mismatch,
             aligner->scoring->gap_open, aligner->scoring->gap_extend);
      printf("\n");
  }
}


