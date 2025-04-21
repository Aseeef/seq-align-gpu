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
#include <limits.h> // INT_MAX
#include <string.h> // memset
#include <ctype.h> // tolower

#include "alignment_scoring.h"
#include "alignment_macros.h"

void scoring_init(scoring_t* scoring,
                  int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool no_gaps_in_a, bool no_gaps_in_b,
                  bool no_mismatches, bool case_sensitive)
{
  // Gap of length 1 has penalty (gap_open+gap_extend)
  // of length N: (gap_open + gap_extend*N)
  scoring->gap_open = gap_open;
  scoring->gap_extend = gap_extend;

  scoring->no_gaps_in_a = no_gaps_in_a;
  scoring->no_gaps_in_b = no_gaps_in_b;
  scoring->no_mismatches = no_mismatches;

  scoring->use_match_mismatch = 1;
  scoring->match = match;
  scoring->mismatch = mismatch;

  scoring->case_sensitive = case_sensitive;

  memset(scoring->wildcards, 0, sizeof(scoring->wildcards));
  memset(scoring->swap_set, 0, sizeof(scoring->swap_set));

  scoring->min_penalty = MIN2(match, mismatch);
  scoring->max_penalty = MAX2(match, mismatch);
  if(!no_gaps_in_a || !no_gaps_in_b) {
    scoring->min_penalty = MIN3(scoring->min_penalty,gap_open+gap_extend,gap_extend);
    scoring->max_penalty = MAX3(scoring->max_penalty,gap_open+gap_extend,gap_extend);
  }
}

void scoring_add_wildcard(scoring_t* scoring, char c, int score)
{
  if(!scoring->case_sensitive) c = tolower(c);
  set_wildcard_bit(scoring,c);
  scoring->wildscores[(size_t)c] = score;
  scoring->min_penalty = MIN2(scoring->min_penalty, score);
  scoring->max_penalty = MAX2(scoring->max_penalty, score);
}

void scoring_add_mutation(scoring_t* scoring, char a, char b, int score)
{
  scoring->swap_scores[(size_t)a][(size_t)b] = score;
  set_swap_bit(scoring,a,b);
  scoring->min_penalty = MIN2(scoring->min_penalty, score);
  scoring->max_penalty = MAX2(scoring->max_penalty, score);
}


// a, b must be lowercase if !scoring->case_sensitive
static char _scoring_check_wildcards(const scoring_t* scoring, char a, char b,
                                     int* score)
{
  // Check if either characters are wildcards
  int tmp_score = INT_MAX;
  if(get_wildcard_bit(scoring,a)) tmp_score = scoring->wildscores[(size_t)a];
  if(get_wildcard_bit(scoring,b)) tmp_score = MIN2(scoring->wildscores[(size_t)b],tmp_score);
  if(tmp_score != INT_MAX) {
    *score = tmp_score;
    return 1;
  }

  *score = 0;
  return 0;
}

// Considered match if lc(a)==lc(b) or if a or b are wildcards
// Always sets score and is_match
void scoring_lookup(const scoring_t* scoring, char a, char b,
                    int *score, bool *is_match)
{
  if(!scoring->case_sensitive)
  {
    a = tolower(a);
    b = tolower(b);
  }

  //#ifdef SEQ_ALIGN_VERBOSE
  //printf(" scoring_lookup(%c,%c)\n", a, b);
  //#endif

  *is_match = (a == b);

  if(scoring->no_mismatches && !*is_match)
  {
    // Check wildcards
    *is_match = _scoring_check_wildcards(scoring, a, b, score);
    return;
  }

  // Look up in table
  if(get_swap_bit(scoring,a,b))
  {
    *score = scoring->swap_scores[(size_t)a][(size_t)b];
    return;
  }

  // Check wildcards
  // Wildcards are used in the order they are given
  // e.g. if we specify '--wildcard X 2 --wildcard Y 3' X:Y align with score 2
  if(_scoring_check_wildcards(scoring, a, b, score))
  {
    *is_match = 1;
    return;
  }

  // Use match/mismatch
  if(scoring->use_match_mismatch)
  {
    *score = (*is_match ? scoring->match : scoring->mismatch);
    return;
  }

  // Error
  fprintf(stderr, "Error: Unknown character pair (%c,%c) and "
                   "match/mismatch have not been set\n", a, b);
  exit(EXIT_FAILURE);
}

// Default
void scoring_system_default(scoring_t* scoring)
{
  int match_default = 1;
  int mismatch_default = -2;
  int gap_open_default = -4;
  int gap_extend_default = -1;

  // case_sensitive = 0
  scoring_init(scoring, match_default, mismatch_default,
               gap_open_default, gap_extend_default,
               0, 0, 0, 0);
}
