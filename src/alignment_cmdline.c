/*
 alignment_cmdline.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// Turn on debugging output by defining SEQ_ALIGN_VERBOSE
//#define SEQ_ALIGN_VERBOSE

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _DEFAULT_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <limits.h> // INT_MIN
#include <stdarg.h> // for va_list
#include <time.h>
#include <omp.h>

#include "seq_file/seq_file.h"

#include "alignment_cmdline.h"

#include "alignment_scoring_load.h"
#include "alignment_scoring.h"

char parse_entire_score_t(char *str, score_t *result) {
    if (sizeof(score_t) == sizeof(int)) {
        return parse_entire_int(str, result);
    } else if (sizeof(score_t) == sizeof(short)) {
        return parse_entire_short(str, result);
    } else {
        // shouldn't happen
        fprintf(stderr, "Error: sizeof(score_t) is not int or short\n");
        exit(EXIT_FAILURE);
    }
}

char parse_entire_short(char *str, short *result) {
    size_t len = strlen(str);

    char *strtol_last_char_ptr = str;
    long tmp = strtol(str, &strtol_last_char_ptr, 10);

    if (tmp > SHRT_MAX || tmp < SHRT_MIN || strtol_last_char_ptr != str + len) {
        return 0;
    } else {
        *result = (short) tmp;
        return 1;
    }
}

char parse_entire_ushort(char *str, unsigned short *result) {
    size_t len = strlen(str);

    char *strtol_last_char_ptr = str;
    long tmp = strtol(str, &strtol_last_char_ptr, 10);

    if (tmp > USHRT_MAX || strtol_last_char_ptr != str + len) {
        return 0;
    } else {
        *result = (short) tmp;
        return 1;
    }
}

char parse_entire_int(char *str, int *result) {
    size_t len = strlen(str);

    char *strtol_last_char_ptr = str;
    long tmp = strtol(str, &strtol_last_char_ptr, 10);

    if (tmp > INT_MAX || tmp < INT_MIN || strtol_last_char_ptr != str + len) {
        return 0;
    } else {
        *result = (int) tmp;
        return 1;
    }
}

char parse_entire_uint(char *str, unsigned int *result) {
    size_t len = strlen(str);

    char *strtol_last_char_ptr = str;
    unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

    if (tmp > UINT_MAX || strtol_last_char_ptr != str + len) {
        return 0;
    } else {
        *result = (unsigned int) tmp;
        return 1;
    }
}

static void print_usage(enum SeqAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...)
__attribute__((format(printf, 4, 5)))
__attribute__((noreturn));

static void print_usage(enum SeqAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...) {
    if (errfmt != NULL) {
        fprintf(stderr, "Error: ");
        va_list argptr;
        va_start(argptr, errfmt);
        vfprintf(stderr, errfmt, argptr);
        va_end(argptr);

        if (errfmt[strlen(errfmt) - 1] != '\n') fprintf(stderr, "\n");
    }

    fprintf(stderr, "usage: %s [OPTIONS] [seq1 seq2]\n", cmdstr);

    fprintf(stderr,
            "  %s optimal %s alignment (maximises score).  \n"
            "  Takes a pair of sequences on the command line, or can read from a\n"
            "  file and from sequence piped in.  Can read gzip files, FASTA and FASTQ.\n\n",
            cmd_type == SEQ_ALIGN_SW_CMD ? "Smith-Waterman" : "Needleman-Wunsch",
            cmd_type == SEQ_ALIGN_SW_CMD ? "local" : "global");

    fprintf(stderr,
            "  OPTIONS:\n"
            "    --file <file>        Sequence file reading with gzip support - read two\n"
            "                         sequences at a time and align them\n"
            "    --files <f1> <f2>    Read one sequence from each file to align at one time\n"
            "    --stdin              Read from STDIN (same as '--file -')\n"
            "\n"
            "    --case_sensitive     Use case sensitive character comparison [default: off]\n"
            "\n"
            "    --match <score>      [default: %i]\n"
            "    --mismatch <score>   [default: %i]\n"
            "    --gapopen <score>    [default: %i]\n"
            "    --gapextend <score>  [default: %i]\n"
            "\n"
            // TODO: worth refractoring to ensure ppl always pass in what i want
            "    --substitution_matrix <file>  see details for formatting\n\n",
            defaults[0], defaults[1],
            defaults[2], defaults[3]);

    if (cmd_type == SEQ_ALIGN_SW_CMD) {
        // SW specific
        fprintf(stderr,
                "    --minscore <score>   Minimum required score\n"
                "                         [default: match * MAX(0.2 * length, 2)]\n"
                "\n"
                "    --printseq           Print sequences before local alignments\n");
    }

    fprintf(stderr,
            "    --printmatrices      Print dynamic programming matrices\n"
            "    --printfasta         Print fasta header lines\n"
            "    --pretty             Print with a descriptor line\n"
            "    --colour             Print with colour\n"
            "\n");

    printf(
        "\n"
        " DETAILS:\n"
        "  * For help choosing scoring, see the README file. \n"
        "  * Gap (of length N) penalty is: (open+N*extend)\n"
        "  * To do alignment without affine gap penalty, set '--gapopen 0'.\n"
        "  * Scoring files should be matrices, with entries separated by a single\n"
        "    character or whitespace. See files in the 'scores' directory for examples.\n"
        "\n"
        "  turner.isaac@gmail.com  (compiled: %s %s)\n", __DATE__, __TIME__);

    exit(EXIT_FAILURE);
}

void cmdline_free(cmdline_t *cmd) {
    free(cmd);
}

#define usage(fmt,...) print_usage(cmd_type,defaults,argv[0],fmt, ##__VA_ARGS__)

cmdline_t *cmdline_new(int argc, char **argv, scoring_t *scoring,
                       enum SeqAlignCmdType cmd_type) {
    cmdline_t *cmd = calloc(1, sizeof(cmdline_t));

    // Store defaults
    score_t defaults[4] = {
        scoring->match, scoring->mismatch,
        scoring->gap_open, scoring->gap_extend
    };

    if (argc == 1)
        usage(NULL);

    // First run through arguments to set up case_sensitive and scoring system

    // case sensitive needs to be dealt with first
    // (is is used to construct hash table for swap_scores)
    char substitutions_set = 0, match_set = 0, mismatch_set = 0;

    int argi;
    for (argi = 1; argi < argc; argi++) {
        if (strcasecmp(argv[argi], "--help") == 0 ||
            strcasecmp(argv[argi], "-help") == 0 ||
            strcasecmp(argv[argi], "-h") == 0) {
            usage(NULL);
        }
    }

    for (argi = 1; argi < argc; argi++) {
        if (argv[argi][0] == '-') {
            if (strcasecmp(argv[argi], "--case_sensitive") == 0) {
                // Already dealt with
                //case_sensitive = true;
            } else if (strcasecmp(argv[argi], "--printseq") == 0) {
                if (cmd_type != SEQ_ALIGN_SW_CMD)
                    usage("--printseq only valid with Smith-Waterman");
                cmd->print_seq = true;
            } else if (strcasecmp(argv[argi], "--printmatrices") == 0) {
                cmd->print_matrices = true;
            } else if (strcasecmp(argv[argi], "--printfasta") == 0) {
                cmd->print_fasta = true;
            } else if (strcasecmp(argv[argi], "--pretty") == 0) {
                cmd->print_pretty = true;
            } else if (strcasecmp(argv[argi], "--colour") == 0) {
                cmd->print_colour = true;
            } else if (strcasecmp(argv[argi], "--stdin") == 0) {
                // Similar to --file argument below
                cmdline_set_files(cmd, "", NULL);
                cmd->interactive = true;
            } else if (argi == argc - 1) {
                // All the remaining options take an extra argument
                usage("Unknown argument without parameter: %s", argv[argi]);
            } else if (strcasecmp(argv[argi], "--scoring") == 0) {
                // This handled above
                argi++; // took an argument
            } else if (strcasecmp(argv[argi], "--substitution_matrix") == 0) {
                gzFile sub_matrix_file = gzopen(argv[argi + 1], "r");
                if (sub_matrix_file == NULL)
                    usage("Couldn't read: %s", argv[argi+1]);

                align_scoring_load_matrix(sub_matrix_file, argv[argi + 1],
                                          scoring, true);

                gzclose(sub_matrix_file);
                substitutions_set = true;

                argi++; // took an argument
            } else if (strcasecmp(argv[argi], "--match") == 0) {
                if (!parse_entire_int(argv[argi + 1], &scoring->match)) {
                    usage("Invalid --match argument ('%s') must be an int", argv[argi+1]);
                }

                match_set = true;
                argi++; // took an argument
            } else if (strcasecmp(argv[argi], "--mismatch") == 0) {
                if (!parse_entire_int(argv[argi + 1], &scoring->mismatch)) {
                    usage("Invalid --mismatch argument ('%s') must be an int", argv[argi+1]);
                }

                mismatch_set = true;
                argi++; // took an argument
            } else if (strcasecmp(argv[argi], "--gapopen") == 0) {
                if (!parse_entire_score_t(argv[argi + 1], &scoring->gap_open)) {
                    usage("Invalid --gapopen argument ('%s') must be an int", argv[argi+1]);
                }

                argi++; // took an argument
            } else if (strcasecmp(argv[argi], "--gapextend") == 0) {
                if (!parse_entire_score_t(argv[argi + 1], &scoring->gap_extend)) {
                    usage("Invalid --gapextend argument ('%s') must be an int",
                          argv[argi+1]);
                }

                argi++; // took an argument
            } else if (strcasecmp(argv[argi], "--file") == 0) {
                cmdline_set_files(cmd, argv[argi + 1], NULL);
                argi++; // took an argument
            }
            // Remaining options take two arguments but check themselves
            else if (strcasecmp(argv[argi], "--files") == 0) {
                printf("Query File=%s and Database File=%s\n", argv[argi + 1], argv[argi + 2]);
                if (argi >= argc - 2) {
                    usage("--files option takes 2 arguments");
                } else if (strcmp(argv[argi + 1], "-") == 0 && strcmp(argv[argi + 2], "-") == 0) {
                    // Read both from stdin
                    cmdline_set_files(cmd, argv[argi + 1], NULL);
                } else {
                    cmdline_set_files(cmd, argv[argi + 1], argv[argi + 2]);
                }

                argi += 2; // took two arguments
            } else
                usage("Unknown argument '%s'", argv[argi]);
        } else {
            if (argc - argi != 2)
                usage("Unknown options: '%s'", argv[argi]);
            break;
        }
    }

    if (substitutions_set && !match_set) {
        // if substitution table set and not match/mismatch
        scoring->use_match_mismatch = 0;
    }

    if (scoring->use_match_mismatch && scoring->match < scoring->mismatch) {
        usage("Match value should not be less than mismatch penalty");
    }

    if (cmd->file_path1 == NULL || cmd->file_path2 == NULL) {
        usage("No input specified");
    }

    return cmd;
}


void cmdline_set_files(cmdline_t *cmd, char *query, char *database) {
    cmd->file_path1 = query;
    cmd->file_path2 = database;
}

char *cmdline_get_file1(cmdline_t *cmd) {
    return cmd->file_path1;
}

char *cmdline_get_file2(cmdline_t *cmd) {
    return cmd->file_path2;
}

static double interval(struct timespec start, struct timespec end) {
    struct timespec temp;
    temp.tv_sec = end.tv_sec - start.tv_sec;
    temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    if (temp.tv_nsec < 0) {
        temp.tv_sec = temp.tv_sec - 1;
        temp.tv_nsec = temp.tv_nsec + 1000000000;
    }
    return (((double) temp.tv_sec) + ((double) temp.tv_nsec) * 1.0e-9);
}

static seq_file_t *open_seq_file(const char *path, bool use_zlib) {
    return (strcmp(path, "-") != 0 || use_zlib)
               ? seq_open(path)
               : seq_dopen(fileno(stdin), false, false, 0);
}

#define BATCH_SIZE_FACTOR 64

void align_from_query_and_db(const char *query_path, const char *db_path, scoring_t *scoring,
                             void (print_alignment)(aligner_t * aligner, size_t total_cnt),
                             bool use_zlib) {
    int num_threads = omp_get_max_threads();
    int max_batch_size = num_threads * BATCH_SIZE_FACTOR;

    size_t VECTOR_SIZE = 32 / sizeof(score_t);
    seq_file_t *query_file, *db_file;
    struct timespec time_start, time_stop;
    size_t i;

    // Open query file
    if ((query_file = open_seq_file(query_path, use_zlib)) == NULL) {
        fprintf(stderr, "Error: couldn't open query file %s\n", query_path);
        fflush(stderr);
        return;
    }

    // Open database file
    if ((db_file = open_seq_file(db_path, use_zlib)) == NULL) {
        fprintf(stderr, "Error: couldn't open database file %s\n", db_path);
        fflush(stderr);
        seq_close(query_file);
        return;
    }

    // Read the single query sequence
    read_t query_read;
    seq_read_alloc(&query_read);

    if (seq_read(query_file, &query_read) <= 0) {
        fprintf(stderr, "Error: Query file %s is empty or invalid\n", query_path);
        fflush(stderr);
        seq_close(query_file);
        seq_close(db_file);
        seq_read_dealloc(&query_read);
        return;
    }

    assert(query_read.name.end != 0);
    char *query_seq = query_read.seq.b;
    char *query_fasta = query_read.name.b;
    size_t query_seq_len = query_read.seq.end;
    assert(query_seq_len > 0);
    // characters are converted into indexes for table lookup
    int8_t *query_indexes = aligned_alloc(32, query_seq_len * sizeof(int8_t));

    // Replace unknown characters in query with an X
    for (i = 0; i < query_seq_len; i++) {
        query_indexes[i] = letters_to_index(query_seq[i]);
        if (!get_swap_bit(scoring, query_indexes[i], query_indexes[i])) {
            query_indexes[i] = letters_to_index('X');
        }
    }

    // Read database sequences and align each with the query
    read_t db_read;
    seq_read_alloc(&db_read);

    // Allocate aligners
    aligner_t **aligners = malloc(sizeof(aligner_t*) * max_batch_size);
    // Init aligners to null
    for (i = 0; i < max_batch_size; i++) {
        aligners[i] = NULL;
    }

    int8_t *db_seq_index_vec_batch = NULL;
    char **db_seq_vec_batch = NULL;
    char **db_fasta_vec_batch = NULL;
    size_t vec_elem_cnt = 0;
    size_t batch_cnt = 0;
    size_t total_cnt = 0;

    bool len_set = false;
    // the max length of the query seq with in the batch
    size_t max_seq_len_in_vec = 0;

    double total_time = 0;

    int read_status = seq_read(db_file, &db_read);
    while (read_status > 0) {
        assert(db_read.name.end != 0);

        char *seq_b = db_read.seq.b;
        size_t seq_b_len = db_read.seq.end;

        if (!len_set) {
            len_set = true;
            // since db is sorted from longest to shortest, first item in batch
            // will always be the longest.
            max_seq_len_in_vec = seq_b_len;
            db_seq_index_vec_batch = aligned_alloc(32, max_seq_len_in_vec * VECTOR_SIZE * sizeof(int8_t));
            db_seq_vec_batch = malloc(sizeof(char *) * VECTOR_SIZE);
            db_fasta_vec_batch = malloc(sizeof(char *) * VECTOR_SIZE);
        }

        assert(max_seq_len_in_vec >= seq_b_len);
        assert(db_seq_index_vec_batch != NULL);
        assert(db_seq_vec_batch != NULL);
        assert(db_fasta_vec_batch != NULL);

        for (i = 0; i < seq_b_len; i++) {
            db_seq_index_vec_batch[i * VECTOR_SIZE + vec_elem_cnt] = letters_to_index(seq_b[i]);
        }
        // characters that are too long are matched with *
        for (i = seq_b_len; i < max_seq_len_in_vec; i++) {
            db_seq_index_vec_batch[i * VECTOR_SIZE + vec_elem_cnt] = letters_to_index('*');
        }
        db_seq_vec_batch[vec_elem_cnt] = strdup(db_read.seq.b);
        db_fasta_vec_batch[vec_elem_cnt] = strdup(db_read.name.b);

        total_cnt++;
        vec_elem_cnt++;

        read_status = seq_read(db_file, &db_read);

        if (vec_elem_cnt == VECTOR_SIZE) {
            assert(query_seq_len != 0);
            assert(max_seq_len_in_vec != 0);
            assert(db_seq_vec_batch != NULL);
            assert(db_fasta_vec_batch != NULL);
            assert(db_seq_index_vec_batch != NULL);


            if (aligners[batch_cnt] == NULL) {
                aligners[batch_cnt] = aligner_create(
                    query_seq,
                    db_seq_vec_batch,
                    query_fasta,
                    db_fasta_vec_batch,
                    query_indexes,
                    db_seq_index_vec_batch,
                    query_seq_len,
                    max_seq_len_in_vec,
                    VECTOR_SIZE,
                    scoring);
            } else {
                aligner_update(
                    aligners[batch_cnt],
                    query_seq,
                    db_seq_vec_batch,
                    query_fasta,
                    db_fasta_vec_batch,
                    query_indexes,
                    db_seq_index_vec_batch,
                    query_seq_len,
                    max_seq_len_in_vec,
                    VECTOR_SIZE,
                    scoring);
            }

            // reset variables for the next vector batch
            len_set = false;
            vec_elem_cnt = 0;
            max_seq_len_in_vec = 0;
            // increment batch
            batch_cnt++;

            if (batch_cnt == max_batch_size || read_status <= 0) {

                clock_gettime(CLOCK_REALTIME, &time_start);
#pragma omp parallel for schedule(static, 1)
                for (i = 0; i < batch_cnt; i++) {
                    alignment_fill_matrices(aligners[i]);
                }
                clock_gettime(CLOCK_REALTIME, &time_stop);
                total_time += interval(time_start, time_stop);

                for (i = 0; i < batch_cnt; i++) {
                    print_alignment(aligners[i], total_cnt - (batch_cnt * VECTOR_SIZE) + (i * VECTOR_SIZE));
                    // cleanup data specific to this batch since it wont be needed again
                    free(aligners[i]->seq_b_batch_indexes);
                    free(aligners[i]->seq_b_str_batch);
                    free(aligners[i]->seq_b_fasta_batch);
                    aligners[i]->seq_b_batch_indexes = NULL;
                    aligners[i]->seq_b_str_batch = NULL;
                    aligners[i]->seq_b_fasta_batch = NULL;
                }

                // reset batch cnt
                batch_cnt = 0;
            }

        }
    }

    printf("Total time: %f\n", total_time);

    // Close files and free memory
    seq_close(query_file);
    seq_close(db_file);
    seq_read_dealloc(&query_read);
    seq_read_dealloc(&db_read);
    free(query_indexes);
    free(aligners);
}
