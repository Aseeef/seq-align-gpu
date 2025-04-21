/*
 examples/sw_example.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: May 2013
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "smith_waterman.h"

void align(char* seq_a, char* seq_b)
{
  // Variables to store alignment result
  sw_aligner_t *sw = smith_waterman_new();

  // Decide on scoring
  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  // Compare character case-sensitively (usually set to 0 for DNA etc)
  char case_sensitive = 0;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend, case_sensitive);

  // Add some special cases
  // x -> y means x in seq1 changing to y in seq2
  scoring_add_mutation(&scoring, 'a', 'c', -2); // a -> c give substitution score -2
  scoring_add_mutation(&scoring, 'c', 'a', -1); // c -> a give substitution score -1

  // We could also prohibit the aligning of characters not given as special cases
  // scoring.use_match_mismatch = 0;

  smith_waterman_align(seq_a, seq_b, &scoring, sw);

  // Free memory for storing alignment results
  smith_waterman_free(sw);
}

int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    printf("usage: ./sw_example <seq1> <seq2>\n");
    exit(EXIT_FAILURE);
  }

  align(argv[1], argv[2]);
  exit(EXIT_SUCCESS);
}
