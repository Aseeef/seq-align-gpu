/*
 smith_waterman.h
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

#ifndef SMITH_WATERMAN_HEADER_SEEN
#define SMITH_WATERMAN_HEADER_SEEN

#include "seq_align.h"
#include "alignment.h"

typedef struct sw_aligner_t sw_aligner_t;

#ifdef __cplusplus
extern "C" {
#endif

sw_aligner_t *smith_waterman_new();
void smith_waterman_free(sw_aligner_t *sw_aligner);

aligner_t* smith_waterman_get_aligner(sw_aligner_t *sw);

/*
 Do not alter seq_a, seq_b or scoring whilst calling this method
 or between calls to smith_waterman_get_hit
*/
void smith_waterman_align_batch(char *seq_a, char **seq_b_batch,
                                score_t * seq_a_indexes, score_t * seq_b_batch_indexes,
                                size_t seq_b_max_batch_len, size_t batch_size,
                                const scoring_t *scoring, sw_aligner_t *sw);

#ifdef __cplusplus
}
#endif

#endif
