/*
 smith_waterman.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// Turn on debugging output by defining SEQ_ALIGN_VERBOSE
//#define SEQ_ALIGN_VERBOSE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sort_r/sort_r.h"

#include "smith_waterman.h"
#include "alignment_macros.h"

// Store alignment here
struct sw_aligner_t
{
  aligner_t aligner;
};

sw_aligner_t* smith_waterman_new()
{
  sw_aligner_t *sw = calloc(1, sizeof(sw_aligner_t));
  return sw;
}

void smith_waterman_free(sw_aligner_t *sw)
{
  aligner_destroy(&(sw->aligner));
  free(sw);
}

aligner_t* smith_waterman_get_aligner(sw_aligner_t *sw)
{
    return &sw->aligner;
}

void smith_waterman_align_batch(const char *a, const char **b_batch, size_t batch_size, const scoring_t *scoring, sw_aligner_t *sw) {
    printf("Aligning batch of size %zu\n", batch_size);
    size_t len_A = strlen(a);
    size_t len_B = strlen(b_batch[0]); // yes strings in all batches must be the same size
    aligner_t *aligner = &sw->aligner;
    aligner_align(aligner, a, b_batch, len_A, len_B, batch_size, scoring);
}
