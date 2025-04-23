/*
 alignment_cmdline.h
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

#ifndef ALIGNMENT_CMDLINE_HEADER_SEEN
#define ALIGNMENT_CMDLINE_HEADER_SEEN

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdarg.h> // required for va_list
#include <stdbool.h>
#include "seq_file/seq_file.h"
#include "alignment.h"

enum SeqAlignCmdType {SEQ_ALIGN_SW_CMD};

typedef struct
{
  // file inputs
  char *file_path1, *file_path2;

  // All values initially 0
  bool case_sensitive;
  int match, mismatch, gap_open, gap_extend;

  // SW specific
  unsigned int max_hits_per_alignment;
  bool max_hits_per_alignment_set;
  bool print_seq;

  // NW specific?
  bool print_matrices;

  // Turns off zlib for stdin
  bool interactive;

  // General output
  bool print_fasta, print_pretty, print_colour;

  // Pair of sequences to align
  const char *seq1, *seq2;
} cmdline_t;

char parse_entire_int(char *str, int *result);
char parse_entire_uint(char *str, unsigned int *result);

cmdline_t* cmdline_new(int argc, char **argv, scoring_t *scoring,
                       enum SeqAlignCmdType cmd_type);
void cmdline_free(cmdline_t* cmd);

void cmdline_set_files(cmdline_t* cmd, char* p1, char* p2);
char* cmdline_get_file1(cmdline_t* cmd);
char* cmdline_get_file2(cmdline_t* cmd);


void align_from_query_and_db(const char *query_path, const char *db_path,
                              void (align)(size_t batch_size, char *query_seq, char *db_seq_batch,
                              int seq_b_len, const char *query_name, const char **db_name),
                              bool use_zlib);

#endif
