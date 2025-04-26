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

#include "smith_waterman.h"

cmdline_t *cmd;
scoring_t scoring;

// Alignment results stored here
sw_aligner_t *sw;

static void sw_set_default_scoring() {
    scoring_system_default(&scoring);

    // Change slightly
    scoring.match = 2;
    scoring.mismatch = -2;
    scoring.gap_open = -2;
    scoring.gap_extend = -1;
}

// Align two sequences against each other to find local alignments between them
void align_batch(size_t batch_size, char * query, char ** db_batch,
                 int32_t * query_indexes, int32_t * db_seq_index_batch, size_t query_len, size_t batch_max_len,
                 const char *seq_a_name, const char **seq_b_name_batch) {

    // Check query has length > 0
    if (query[0] == '\0') {
        fprintf(stderr, "Error: The query must have length > 0\n");
        fflush(stderr);

        if (cmd->print_fasta && seq_a_name != NULL) {
            fprintf(stderr, "%s\n", seq_a_name);
        }

        fflush(stderr);

        return;
    }

    smith_waterman_align_batch(query, db_batch,
                               query_indexes, db_seq_index_batch,
                               query_len, batch_max_len, batch_size,
                               &scoring, sw);

    aligner_t *aligner = smith_waterman_get_aligner(sw);

    if (cmd->print_matrices) {
        alignment_print_matrices(aligner, batch_size);
    }

    // seqA
    if (cmd->print_fasta && seq_a_name != NULL) {
        fputs(seq_a_name, stdout);
        putc('\n', stdout);
    }

    if (cmd->print_seq) {
        fputs(query, stdout);
        putc('\n', stdout);
    }

    printf("------\n");

    for (size_t b = 0; b < batch_size; b++) {
        // seqB
        if (cmd->print_fasta && seq_b_name_batch[b] != NULL) {
            fputs(seq_b_name_batch[b], stdout);
            putc('\n', stdout);
        }

        if (cmd->print_seq) {
            fputs(db_batch[b], stdout);
            putc('\n', stdout);
        }

        putc('\n', stdout);

        fflush(stdout);

        printf("score: %i\n\n", aligner->max_scores[b]);
    }

    fputs("==\n", stdout);
    fflush(stdout);
}

int main(int argc, char *argv[]) {
#ifdef SEQ_ALIGN_VERBOSE
    printf("VERBOSE: on\n");
#endif

    sw_set_default_scoring();
    cmd = cmdline_new(argc, argv, &scoring, SEQ_ALIGN_SW_CMD);

    // Align!
    sw = smith_waterman_new();

    // Align from files
    const char *query_file = cmdline_get_file1(cmd);
    const char *db_file = cmdline_get_file2(cmd);

    if (query_file != NULL && db_file != NULL) {
        align_from_query_and_db(query_file, db_file, &scoring, &align_batch, !cmd->interactive);
    } else {
        fprintf(stderr, "Error: Both query and database files must be provided\n");
        fflush(stderr);
    }

    // Free memory for storing alignment results
    smith_waterman_free(sw);

    cmdline_free(cmd);

    return EXIT_SUCCESS;
}
