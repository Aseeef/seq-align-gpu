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
inline static __m256i scoring_lookup(const scoring_t *scoring, int8_t a_index, int8_t *b_indexes) {
    // base address (the row) of the swap_scores for a_index
    const int8_t *swap_scores = scoring->swap_scores[a_index];
    alignas(32) int16_t indexes[16];
    for (int32_t i = 0; i < 16; i++) {
        indexes[i] = (int16_t) swap_scores[b_indexes[i]];
    }
    return _mm256_load_si256((__m256i *) indexes);
}

const size_t FULL_VECTOR_SIZE = 32 / sizeof(score_t);

// Fill in traceback matrix for an ENTIRE BATCH
void alignment_fill_matrices(aligner_t *aligner) {
    score_t *curr_match_scores = aligner->curr_match_scores;
    score_t *curr_gap_a_scores = aligner->curr_gap_a_scores;
    score_t *curr_gap_b_scores = aligner->curr_gap_b_scores;
    const scoring_t *scoring = aligner->scoring;
    size_t score_width = aligner->score_width;
    size_t score_height = aligner->score_height;
    size_t i, j;

    __m256i gap_open_penalty = _mm256_set1_epi16(scoring->gap_extend + scoring->gap_open);
    __m256i gap_extend_penalty = _mm256_set1_epi16(scoring->gap_extend);
    __m256i min_v = _mm256_setzero_si256();

    // null checks
    assert(curr_match_scores != NULL);
    assert(curr_gap_a_scores != NULL);
    assert(curr_gap_b_scores != NULL);
    assert(scoring != NULL);

    size_t seq_i, seq_j, len_i = score_width - 1, len_j = score_height - 1;
    size_t index, index_right;

    // reset match scores
    __m256i max_scores_vec = _mm256_setzero_si256();

    // reset match and gap matrices
    // The shape is height x width x b
    // I need to reset first col and first row of each batch
    for (i = 0; i < score_width; i++) {
        size_t offset = i * FULL_VECTOR_SIZE;
        _mm256_store_si256((__m256i *) (curr_match_scores + offset), min_v);
        _mm256_store_si256((__m256i *) (curr_gap_b_scores + offset), min_v);
    }
    for (j = 0; j < score_width; j++) {
        size_t offset = j * FULL_VECTOR_SIZE;
        _mm256_store_si256((__m256i *) (curr_gap_a_scores + offset), min_v);
    }


    for (seq_j = 0; seq_j < len_j; seq_j++) {

        // init these to zeros since we know the the left boundary is all zeros
        __m256i match_score_left = _mm256_setzero_si256();
        __m256i gap_a_score_left = _mm256_setzero_si256();
        __m256i gap_b_score_left = _mm256_setzero_si256();

        __m256i match_score_up_left = _mm256_setzero_si256();
        __m256i gap_a_score_up_left = _mm256_setzero_si256();
        __m256i gap_b_score_up_left = _mm256_setzero_si256();

        // Indices (relative to the single row buffer)
        index = FULL_VECTOR_SIZE; // Start calculating column 1
        index_right = (2 * FULL_VECTOR_SIZE);

        for (seq_i = 0; seq_i < len_i; seq_i++) {

            // Make sure to keep the entire table fetched at all times by prefetching
            // everytime the inner loop finishes
            char const* prefetch_addr = (char const*)scoring->swap_scores;
            _mm_prefetch(prefetch_addr, _MM_HINT_T0);
            _mm_prefetch(prefetch_addr + 64, _MM_HINT_T0);
            _mm_prefetch(prefetch_addr + 128, _MM_HINT_T0);
            _mm_prefetch(prefetch_addr + 192, _MM_HINT_T0);

            // substitution penalty
            __m256i substitution_penalty = scoring_lookup(scoring, aligner->seq_a_indexes[seq_i],
                                                          aligner->seq_b_batch_indexes + (seq_j * FULL_VECTOR_SIZE));

            // Prefetch gap_a_scores, gap_b_scores, match_scores at the upcoming index
            char const *gap_a_scores_ptr = (char const *) (curr_gap_a_scores + index_right);
            char const *gap_b_scores_ptr = (char const *) (curr_gap_b_scores + index_right);
            char const *curr_match_scores_ptr = (char const *) (curr_match_scores + index_right);
            _mm_prefetch(gap_a_scores_ptr, _MM_HINT_T0);
            _mm_prefetch(gap_b_scores_ptr, _MM_HINT_T0);
            _mm_prefetch(curr_match_scores_ptr, _MM_HINT_T0);

            // Currently index has the values of the table from the previous iteration of seq_j (i.e. the row)
            // so we gotta cache em before its overwritten because we need this
            __m256i match_score_up = _mm256_load_si256((__m256i *) (curr_match_scores + index));
            __m256i gap_a_score_up = _mm256_load_si256((__m256i *) (curr_gap_a_scores + index));
            __m256i gap_b_score_up = _mm256_load_si256((__m256i *) (curr_gap_b_scores + index));

            // Update match_scores[i][j]
            //          score_t match_score = MAX4(match_scores[index_upleft] + substitution_penalty,
            //                                     gap_a_scores[index_upleft] + substitution_penalty,
            //                                     gap_b_scores[index_upleft] + substitution_penalty,
            //                                     min);
            // H[i][j] = MAX(0, H[i-1][j-1] + substitution_penalty, F[i-1][j-1] + substitution_penalty, E[i-1][j-1] + substitution_penalty)

            __m256i match_score_curr = _mm256_add_epi16(match_score_up_left, substitution_penalty);
            __m256i gap_a_score_val = _mm256_add_epi16(gap_a_score_up_left, substitution_penalty);
            __m256i gap_b_score_val = _mm256_add_epi16(gap_b_score_up_left, substitution_penalty);
            match_score_curr = _mm256_max_epi16(match_score_curr, gap_a_score_val);
            match_score_curr = _mm256_max_epi16(match_score_curr, gap_b_score_val);
            match_score_curr = _mm256_max_epi16(match_score_curr, min_v);

            // update best score
            // equal to: max_scores_vec[i] = (match_score[i] > max_scores_vec[i]) ? match_score[i] : max_scores_vec[i];
            __m256i mask = _mm256_cmpgt_epi16(match_score_curr, max_scores_vec);
            max_scores_vec = _mm256_blendv_epi8(max_scores_vec, match_score_curr, mask);

            // Update gap_a_scores[i][j]
            //          gap_a_scores[index]
            //                  = MAX4(match_scores[index_up] + gap_open_penalty,
            //                         gap_a_scores[index_up] + gap_extend_penalty,
            //                         gap_b_scores[index_up] + gap_open_penalty,
            //                         min);
            // E[i][j] = MAX( 0, H[i-1][j] + gap_open_penalty, E[i-1][j] + gap_extend_penalty , F[i-1][j] + gap_open_penalty )
            __m256i match_score_val = _mm256_add_epi16(match_score_up, gap_open_penalty);
            gap_a_score_val = _mm256_add_epi16(gap_a_score_up, gap_extend_penalty);
            gap_b_score_val = _mm256_add_epi16(gap_b_score_up, gap_open_penalty);
            __m256i gap_a_score_curr = _mm256_max_epi16(match_score_val, gap_a_score_val);
            gap_a_score_curr = _mm256_max_epi16(gap_a_score_curr, gap_b_score_val);
            gap_a_score_curr = _mm256_max_epi16(gap_a_score_curr, min_v);

            // Update gap_b_scores[i][j]
            //          gap_b_scores[index]
            //                  = MAX4(match_scores[index_left] + gap_open_penalty,
            //                         gap_a_scores[index_left] + gap_open_penalty,
            //                         gap_b_scores[index_left] + gap_extend_penalty,
            //                         min);
            // F[i][j] = MAX( 0, H[i][j-1] + gap_open_penalty, E[i][j-1] + gap_open_penalty , F[i][j-1] + gap_extend_penalty )
            match_score_val = _mm256_add_epi16(match_score_left, gap_open_penalty);
            gap_a_score_val = _mm256_add_epi16(gap_a_score_left, gap_open_penalty);
            gap_b_score_val = _mm256_add_epi16(gap_b_score_left, gap_extend_penalty);
            __m256i gap_b_score_curr = _mm256_max_epi16(match_score_val, gap_a_score_val);
            gap_b_score_curr = _mm256_max_epi16(gap_b_score_curr, gap_b_score_val);
            gap_b_score_curr = _mm256_max_epi16(gap_b_score_curr, min_v);

            // Update the buffers
            _mm256_store_si256((__m256i *) (curr_match_scores + index), match_score_curr);
            _mm256_store_si256((__m256i *) (curr_gap_a_scores + index), gap_a_score_curr);
            _mm256_store_si256((__m256i *) (curr_gap_b_scores + index), gap_b_score_curr);


            match_score_up_left = match_score_up;
            gap_a_score_up_left = gap_a_score_up;
            gap_b_score_up_left = gap_b_score_up;


            match_score_left = match_score_curr;
            gap_a_score_left = gap_a_score_curr;
            gap_b_score_left = gap_b_score_curr;

            // inc indexes
            index += FULL_VECTOR_SIZE;
            index_right += FULL_VECTOR_SIZE;
        }
    }

    // put back the max scores in this batch
    assert(aligner->max_scores != NULL);
    _mm256_storeu_si256((__m256i *) (aligner->max_scores), max_scores_vec);
}

// Note: len_b must be same for all batches
void aligner_update(aligner_t *aligner,
                    char *seq_a_str, char **seq_b_str_batch,
                    char *seq_a_fasta, char **seq_b_fasta_batch,
                    int8_t *seq_a_indexes, int8_t *seq_b_batch_indexes,
                    size_t len_a, size_t len_b, size_t vector_size,
                    const scoring_t *scoring) {
    aligner->scoring = scoring;
    aligner->seq_a_str = seq_a_str;
    aligner->seq_b_str_batch = seq_b_str_batch;
    aligner->seq_a_indexes = seq_a_indexes;
    aligner->seq_b_batch_indexes = seq_b_batch_indexes;
    aligner->seq_a_fasta = seq_a_fasta;
    aligner->seq_b_fasta_batch = seq_b_fasta_batch;
    aligner->vector_size = vector_size;
    aligner->score_width = len_a + 1; // for col of all zeros
    aligner->score_height = len_b + 1; // for the row of all zeros
}

aligner_t *aligner_create(char *seq_a_str, char **seq_b_str_batch,
                          char *seq_a_fasta, char **seq_b_fasta_batch,
                          int8_t *seq_a_indexes, int8_t *seq_b_batch_indexes,
                          size_t len_a, size_t len_b, size_t vector_size,
                          const scoring_t *scoring) {
    aligner_t *aligner = malloc(sizeof(aligner_t));
    aligner->scoring = scoring;
    aligner->seq_a_str = seq_a_str;
    aligner->seq_b_str_batch = seq_b_str_batch;
    aligner->seq_a_fasta = seq_a_fasta;
    aligner->seq_b_fasta_batch = seq_b_fasta_batch;
    aligner->seq_a_indexes = seq_a_indexes;
    aligner->seq_b_batch_indexes = seq_b_batch_indexes;
    aligner->vector_size = vector_size;
    aligner->score_width = len_a + 1; // for col of all zeros
    aligner->score_height = len_b + 1; // for the row of all zeros

    aligner->max_scores = aligned_alloc(32, sizeof(score_t) * vector_size);
    // arrays are traversed row by row so h_mem makes sense
    size_t h_mem_size = sizeof(score_t) * aligner->score_width * vector_size;
    aligner->curr_match_scores = aligned_alloc(32, h_mem_size);
    aligner->curr_gap_a_scores = aligned_alloc(32, h_mem_size);
    aligner->curr_gap_b_scores = aligned_alloc(32, h_mem_size);

    return aligner;
}

void aligner_destroy(aligner_t *aligner) {
    free(aligner->curr_match_scores);
    free(aligner->curr_gap_a_scores);
    free(aligner->curr_gap_b_scores);
}


// void alignment_print_matrices(const aligner_t *aligner) {
//     const score_t *match_scores = aligner->match_scores;
//     const score_t *gap_a_scores = aligner->gap_a_scores;
//     const score_t *gap_b_scores = aligner->gap_b_scores;
//
//     size_t i, j, b, vector_size;
//     vector_size = aligner->vector_size;
//
//
//     // Note: Yes, the indexing here is not cache friendly but it's fine, this code is only
//     // used to debug. So I rather not refactor it...
//     for (b = 0; b < vector_size; b++) {
//         printf("batch no: %zi/%zi,\n seq_a: %.*s,\nseq_b: %.*s\n",
//                b + 1, vector_size,
//                (int) aligner->score_width - 1, aligner->seq_a_str,
//                (int) strlen(aligner->seq_b_str_batch[b]), aligner->seq_b_str_batch[b]);
//
//         printf("\n");
//
//         printf("match_scores:\n");
//         for (j = 0; j < aligner->score_height; j++) {
//             printf("%3i:", (int) j);
//             for (i = 0; i < aligner->score_width; i++) {
//                 printf("\t%3i", (int) match_scores[ARR_3D_index(aligner->score_width, FULL_VECTOR_SIZE, j, i, b)]);
//                 //printf("\t%3i", (int)ARR_LOOKUP(match_scores, aligner->score_width, i, j));
//             }
//             putc('\n', stdout);
//         }
//         printf("gap_a_scores:\n");
//         for (j = 0; j < aligner->score_height; j++) {
//             printf("%3i:", (int) j);
//             for (i = 0; i < aligner->score_width; i++) {
//                 printf("\t%3i", (int) gap_a_scores[ARR_3D_index(aligner->score_width, FULL_VECTOR_SIZE, j, i, b)]);
//             }
//             putc('\n', stdout);
//         }
//         printf("gap_b_scores:\n");
//         for (j = 0; j < aligner->score_height; j++) {
//             printf("%3i:", (int) j);
//             for (i = 0; i < aligner->score_width; i++) {
//                 printf("\t%3i", (int) gap_b_scores[ARR_3D_index(aligner->score_width, FULL_VECTOR_SIZE, j, i, b)]);
//             }
//             putc('\n', stdout);
//         }
//
//         printf("match: %i mismatch: %i gapopen: %i gapexend: %i\n",
//                aligner->scoring->match, aligner->scoring->mismatch,
//                aligner->scoring->gap_open, aligner->scoring->gap_extend);
//         printf("\n");
//     }
// }
