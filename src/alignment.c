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

    // null checks
    assert(match_scores != NULL);
    assert(gap_a_scores != NULL);
    assert(gap_b_scores != NULL);
    assert(scoring != NULL);

    size_t seq_i, seq_j, len_i = score_width - 1, len_j = score_height - 1;
    size_t index, index_left, index_right, index_up, index_upleft;

    // reset match scores
    __m256i max_scores_vec = _mm256_setzero_si256();

    // reset match and gap matrices
    // The shape is height x width x b
    // I need to reset first col and first row of each batch
    __m256i zero_v = _mm256_setzero_si256();
    for (i = 0; i < score_width; i++) {
        size_t offset = i * FULL_VECTOR_SIZE;
        _mm256_store_si256((__m256i *) (match_scores + offset), zero_v);
        _mm256_store_si256((__m256i *) (gap_a_scores + offset), zero_v);
        _mm256_store_si256((__m256i *) (gap_b_scores + offset), zero_v);
    }
    for (j = 0; j < score_height; j++) {
        size_t offset = j * score_width * FULL_VECTOR_SIZE;
        _mm256_store_si256((__m256i *) (match_scores + offset), zero_v);
        _mm256_store_si256((__m256i *) (gap_a_scores + offset), zero_v);
        _mm256_store_si256((__m256i *) (gap_b_scores + offset), zero_v);
    }

    // start at position [1][1]
    index_upleft = 0;
    index_up = FULL_VECTOR_SIZE;
    index_left = score_width * FULL_VECTOR_SIZE;
    index = (score_width * FULL_VECTOR_SIZE) + FULL_VECTOR_SIZE;
    index_right = (score_width * FULL_VECTOR_SIZE) + (2 * FULL_VECTOR_SIZE);

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

            char const* prefetch_addr = (char const*)scoring->swap_scores;
            _mm_prefetch(prefetch_addr, _MM_HINT_T0);
            _mm_prefetch(prefetch_addr + 64, _MM_HINT_T0);
            _mm_prefetch(prefetch_addr + 128, _MM_HINT_T0);
            _mm_prefetch(prefetch_addr + 192, _MM_HINT_T0);

            // substitution penalty
            __m256i substitution_penalty = scoring_lookup(scoring, aligner->seq_a_indexes[seq_i],
                                                          aligner->seq_b_batch_indexes + (seq_j * FULL_VECTOR_SIZE));

            // Prefetch gap_a_scores, gap_b_scores, match_scores at the upcoming index
            char const *gap_a_scores_ptr = (char const *) (gap_a_scores + index_right);
            char const *gap_b_scores_ptr = (char const *) (gap_b_scores + index_right);
            char const *match_scores_ptr = (char const *) (match_scores + index_right);
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

            __m256i match_score = _mm256_load_si256((__m256i *) (match_scores + index_upleft));
            match_score = _mm256_add_epi16(match_score, substitution_penalty);

            __m256i gap_a_score = _mm256_load_si256((__m256i *) (gap_a_scores + index_upleft));
            gap_a_score = _mm256_add_epi16(gap_a_score, substitution_penalty);
            __m256i gap_b_score = _mm256_load_si256((__m256i *) (gap_b_scores + index_upleft));
            gap_b_score = _mm256_add_epi16(gap_b_score, substitution_penalty);

            __m256i max_ab = _mm256_max_epi16(gap_a_score, gap_b_score);
            __m256i max_mab = _mm256_max_epi16(match_score, max_ab);
            match_score = _mm256_max_epi16(max_mab, min_v);

            // save
            _mm256_store_si256((__m256i *) (match_scores + index), match_score);

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
            match_score = _mm256_load_si256((__m256i *) (match_scores + index_up));
            match_score = _mm256_add_epi16(match_score, gap_open_penalty);
            gap_a_score = _mm256_load_si256((__m256i *) (gap_a_scores + index_up));
            gap_a_score = _mm256_add_epi16(gap_a_score, gap_extend_penalty);
            gap_b_score = _mm256_load_si256((__m256i *) (gap_b_scores + index_up));
            gap_b_score = _mm256_add_epi16(gap_b_score, gap_open_penalty);
            gap_a_score = _mm256_max_epi16(gap_a_score, match_score);
            gap_a_score = _mm256_max_epi16(gap_a_score, gap_b_score);
            gap_a_score = _mm256_max_epi16(gap_a_score, min_v);
            _mm256_store_si256((__m256i *) (gap_a_scores + index), gap_a_score);

            // Update gap_b_scores[i][j]
            //          gap_b_scores[index]
            //                  = MAX4(match_scores[index_left] + gap_open_penalty,
            //                         gap_a_scores[index_left] + gap_open_penalty,
            //                         gap_b_scores[index_left] + gap_extend_penalty,
            //                         min);
            match_score = _mm256_load_si256((__m256i *) (match_scores + index_left));
            match_score = _mm256_add_epi16(match_score, gap_open_penalty);
            gap_a_score = _mm256_load_si256((__m256i *) (gap_a_scores + index_left));
            gap_a_score = _mm256_add_epi16(gap_a_score, gap_open_penalty);
            gap_b_score = _mm256_load_si256((__m256i *) (gap_b_scores + index_left));
            gap_b_score = _mm256_add_epi16(gap_b_score, gap_extend_penalty);
            gap_b_score = _mm256_max_epi16(gap_b_score, match_score);
            gap_b_score = _mm256_max_epi16(gap_b_score, gap_a_score);
            gap_b_score = _mm256_max_epi16(gap_b_score, min_v);
            _mm256_store_si256((__m256i *) (gap_b_scores + index), gap_b_score);

            index += FULL_VECTOR_SIZE;
            index_left += FULL_VECTOR_SIZE;
            index_up += FULL_VECTOR_SIZE;
            index_upleft += FULL_VECTOR_SIZE;
            index_right += FULL_VECTOR_SIZE;
        }
        // need to increment again to make everyone
        // get to the next row
        index += FULL_VECTOR_SIZE;
        index_left += FULL_VECTOR_SIZE;
        index_up += FULL_VECTOR_SIZE;
        index_upleft += FULL_VECTOR_SIZE;
        index_right += FULL_VECTOR_SIZE;
    }

    // put back the max scores in this batch
    assert(aligner->match_scores != NULL);
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

    // The first allocation is expected to be the largest width*height.
    assert(aligner->capacity == 0 || ROUNDUP2POW(aligner->score_width * aligner->score_height) <= aligner->capacity);
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

    size_t capacity = ROUNDUP2POW(aligner->score_width * aligner->score_height);
    aligner->capacity = capacity;
    aligner->max_scores = aligned_alloc(32, sizeof(score_t) * vector_size);
    size_t mem_size = sizeof(score_t) * capacity * vector_size;
    aligner->match_scores = aligned_alloc(32, mem_size);
    aligner->gap_a_scores = aligned_alloc(32, mem_size);
    aligner->gap_b_scores = aligned_alloc(32, mem_size);

    return aligner;
}

void aligner_destroy(aligner_t *aligner) {
    if (aligner->capacity > 0) {
        free(aligner->match_scores);
        free(aligner->gap_a_scores);
        free(aligner->gap_b_scores);
    }
}


void alignment_print_matrices(const aligner_t *aligner) {
    const score_t *match_scores = aligner->match_scores;
    const score_t *gap_a_scores = aligner->gap_a_scores;
    const score_t *gap_b_scores = aligner->gap_b_scores;

    size_t i, j, b, vector_size;
    vector_size = aligner->vector_size;


    // Note: Yes, the indexing here is not cache friendly but it's fine, this code is only
    // used to debug. So I rather not refactor it...
    for (b = 0; b < vector_size; b++) {
        printf("batch no: %zi/%zi,\n seq_a: %.*s,\nseq_b: %.*s\n",
               b + 1, vector_size,
               (int) aligner->score_width - 1, aligner->seq_a_str,
               (int) strlen(aligner->seq_b_str_batch[b]), aligner->seq_b_str_batch[b]);

        printf("\n");

        printf("match_scores:\n");
        for (j = 0; j < aligner->score_height; j++) {
            printf("%3i:", (int) j);
            for (i = 0; i < aligner->score_width; i++) {
                printf("\t%3i", (int) match_scores[ARR_3D_index(aligner->score_width, FULL_VECTOR_SIZE, j, i, b)]);
                //printf("\t%3i", (int)ARR_LOOKUP(match_scores, aligner->score_width, i, j));
            }
            putc('\n', stdout);
        }
        printf("gap_a_scores:\n");
        for (j = 0; j < aligner->score_height; j++) {
            printf("%3i:", (int) j);
            for (i = 0; i < aligner->score_width; i++) {
                printf("\t%3i", (int) gap_a_scores[ARR_3D_index(aligner->score_width, FULL_VECTOR_SIZE, j, i, b)]);
            }
            putc('\n', stdout);
        }
        printf("gap_b_scores:\n");
        for (j = 0; j < aligner->score_height; j++) {
            printf("%3i:", (int) j);
            for (i = 0; i < aligner->score_width; i++) {
                printf("\t%3i", (int) gap_b_scores[ARR_3D_index(aligner->score_width, FULL_VECTOR_SIZE, j, i, b)]);
            }
            putc('\n', stdout);
        }

        printf("match: %i mismatch: %i gapopen: %i gapexend: %i\n",
               aligner->scoring->match, aligner->scoring->mismatch,
               aligner->scoring->gap_open, aligner->scoring->gap_extend);
        printf("\n");
    }
}
