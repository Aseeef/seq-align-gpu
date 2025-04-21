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

/**
 * Initializes a scoring_t object with the given scoring parameters.
 *
 * @param scoring          Pointer to the scoring_t structure to be initialized.
 * @param match            Score for a match between characters.
 * @param mismatch         Penalty score for a mismatch between characters.
 * @param gap_open         Penalty score for opening a gap.
 * @param gap_extend       Penalty score for extending a gap.
 * @param case_sensitive   If true, scoring is case-sensitive (uppercase and lowercase treated differently).
 */
void scoring_init(scoring_t *scoring,
                  int match, int mismatch,
                  int gap_open, int gap_extend,
                  bool case_sensitive) {
    // Gap of length 1 has penalty (gap_open+gap_extend)
    // of length N: (gap_open + gap_extend*N)
    scoring->gap_open = gap_open;
    scoring->gap_extend = gap_extend;

    scoring->use_match_mismatch = 1;
    scoring->match = match;
    scoring->mismatch = mismatch;

    scoring->case_sensitive = case_sensitive;

    memset(scoring->swap_set, 0, sizeof(scoring->swap_set));

    scoring->min_penalty = MIN2(match, mismatch);
    scoring->max_penalty = MAX2(match, mismatch);
}

/**
 * Adds a mutation score between two specific characters.
 *
 * @param scoring          Pointer to the scoring_t structure.
 * @param a                First character involved in the mutation.
 * @param b                Second character involved in the mutation.
 * @param score            Score for the mutation (alignment between a and b).
 */
void scoring_add_mutation(scoring_t *scoring, char a, char b, int score) {
    scoring->swap_scores[(size_t) a][(size_t) b] = score;
    set_swap_bit(scoring, a, b);
    scoring->min_penalty = MIN2(scoring->min_penalty, score);
    scoring->max_penalty = MAX2(scoring->max_penalty, score);
}

/**
 * Checks if characters a or b are wildcards and fetches the associated score.
 *
 * @param scoring          Pointer to the scoring_t structure.
 * @param a                First character to check for wildcard matching.
 * @param b                Second character to check for wildcard matching.
 * @param score            Pointer to store the score if wildcards match.
 * 
 * @return                 Returns 1 if either character is a wildcard, else 0.
 */
static char _scoring_check_wildcards(const scoring_t *scoring, char a, char b,
                                     int *score) {
    return 0;
}

/**
 * Looks up the score for aligning characters a and b and determines if they match.
 *
 * @param scoring          Pointer to the scoring_t structure.
 * @param a                First character in the alignment.
 * @param b                Second character in the alignment.
 * @param score            Pointer to store the score for aligning a and b.
 * @param is_match         Pointer to store whether a and b are considered a match.
 */
void scoring_lookup(const scoring_t *scoring, char a, char b,
                    int *score, bool *is_match) {
    if (!scoring->case_sensitive) {
        a = tolower(a);
        b = tolower(b);
    }

    //#ifdef SEQ_ALIGN_VERBOSE
    //printf(" scoring_lookup(%c,%c)\n", a, b);
    //#endif

    *is_match = (a == b);

    // Look up in table
    if (get_swap_bit(scoring, a, b)) {
        *score = scoring->swap_scores[(size_t) a][(size_t) b];
        return;
    }

    // Check wildcards
    // Wildcards are used in the order they are given
    // e.g. if we specify '--wildcard X 2 --wildcard Y 3' X:Y align with score 2
    if (_scoring_check_wildcards(scoring, a, b, score)) {
        *is_match = 1;
        return;
    }

    // Use match/mismatch
    if (scoring->use_match_mismatch) {
        *score = (*is_match ? scoring->match : scoring->mismatch);
        return;
    }

    // Error
    fprintf(stderr, "Error: Unknown character pair (%c,%c) and "
                    "match/mismatch have not been set\n", a, b);
    exit(EXIT_FAILURE);
}

/**
 * Sets the scoring_t structure with default scoring parameters.
 *
 * @param scoring          Pointer to the scoring_t structure to be initialized with defaults.
 */
void scoring_system_default(scoring_t *scoring) {
    int match_default = 1;
    int mismatch_default = -2;
    int gap_open_default = -4;
    int gap_extend_default = -1;

    // case_sensitive = 0
    scoring_init(scoring, match_default, mismatch_default,
                 gap_open_default, gap_extend_default,
                 0);
}
