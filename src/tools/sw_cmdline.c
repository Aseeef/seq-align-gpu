/*
 tools/sw_cmdline.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// my utility functions
#include "seq_file/seq_file.h"

// Alignment scoring and loading
#include "alignment_scoring_load.h"
#include "alignment_cmdline.h"
#include "alignment_macros.h"

cmdline_t *cmd;

static void sw_set_default_scoring(scoring_t * scoring) {
    scoring_system_default(scoring);

    // Change slightly
    scoring->match = 2;
    scoring->mismatch = -2;
    scoring->gap_open = -2;
    scoring->gap_extend = -1;
}

// Align two sequences against each other to find local alignments between them
void align_batch(aligner_t * aligner, size_t total_cnt) {

    // Check query has length > 0
    if (aligner->seq_a_str[0] == '\0') {
        fprintf(stderr, "Error: The query must have length > 0\n");
        fflush(stderr);

        if (cmd->print_fasta && aligner->seq_a_fasta != NULL) {
            fprintf(stderr, "%s\n", aligner->seq_a_fasta);
        }

        fflush(stderr);

        return;
    }

    alignment_fill_matrices(aligner);

    if (cmd->print_matrices) {
        //alignment_print_matrices(aligner);
    }

    // seqA
    if (cmd->print_fasta && aligner->seq_a_fasta != NULL) {
        fputs(aligner->seq_a_fasta, stdout);
        putc('\n', stdout);
    }

    if (cmd->print_seq) {
        fputs(aligner->seq_a_str, stdout);
        putc('\n', stdout);
    }

    for (size_t b = 0; b < aligner->vector_size; b++) {
        printf("Entry #%lu:\n", total_cnt + b);
        assert(aligner->seq_b_fasta_batch != NULL);
        // seqB
        if (cmd->print_fasta && aligner->seq_b_fasta_batch[b] != NULL) {
            fputs(aligner->seq_b_fasta_batch[b], stdout);
            putc('\n', stdout);
        }

        if (cmd->print_seq) {
            fputs(aligner->seq_b_str_batch[b], stdout);
            putc('\n', stdout);
        }

        fflush(stdout);

        printf("score: %i\n\n", aligner->max_scores[b]);
    }

    fflush(stdout);
}

int main(int argc, char *argv[]) {
#ifdef SEQ_ALIGN_VERBOSE
    printf("VERBOSE: on\n");
#endif

    scoring_t scoring;
    sw_set_default_scoring(&scoring);
    cmd = cmdline_new(argc, argv, &scoring, SEQ_ALIGN_SW_CMD);

    // Align from files
    const char *query_file = cmdline_get_file1(cmd);
    const char *db_file = cmdline_get_file2(cmd);

    if (query_file != NULL && db_file != NULL) {
        align_from_query_and_db(query_file, db_file, &scoring, &align_batch, !cmd->interactive);
    } else {
        fprintf(stderr, "Error: Both query and database files must be provided\n");
        fflush(stderr);
    }

    cmdline_free(cmd);

    return EXIT_SUCCESS;
}
