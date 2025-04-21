/*
 tools/sw_cmdline.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

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

size_t alignment_index = 0;
bool wait_on_keystroke = 0;

static void sw_set_default_scoring()
{
  scoring_system_default(&scoring);

  // Change slightly
  scoring.match = 2;
  scoring.mismatch = -2;
  scoring.gap_open = -2;
  scoring.gap_extend = -1;
}

// Print one line of an alignment
void print_alignment_part(const char* seq1, const char* seq2,
                          size_t pos, size_t len)
{
  size_t i;
  printf("  ");

  if(cmd->print_colour)
    alignment_colour_print_against(seq1, seq2, scoring.case_sensitive);
  else
    fputs(seq1, stdout);

  printf("  [pos: %li; len: %lu]\n", pos, len);
}

static char get_next_hit()
{
  if(!wait_on_keystroke)
    return 1;

  int r = 0;

  char response = 0;
  char next_hit = 0;

  while(!response)
  {
    printf("next [h]it or [a]lignment: ");
    fflush(stdout);

    while((r = getc(stdin)) != -1 && r != '\n' && r != '\r')
    {
      if(r == 'h' || r == 'H')
      {
        next_hit = 1;
        response = 1;
      }
      else if(r == 'a' || r == 'A')
      {
        next_hit = 0;
        response = 1;
      }
    }

    if(r == -1)
    {
      // We're done here
      putc('\n', stdout);
      exit(EXIT_SUCCESS);
    }
  }

  return next_hit;
}

// Align two sequences against each other to find local alignments between them
void align(const char *seq_a, const char *seq_b,
           const char *seq_a_name, const char *seq_b_name)
{
  if((seq_a_name != NULL || seq_b_name != NULL) && wait_on_keystroke)
  {
    fprintf(stderr, "Error: Interactive input takes seq only "
                    "(no FASTA/FASTQ) '%s:%s'\n", seq_a_name, seq_b_name);
    fflush(stderr);
    exit(EXIT_FAILURE);
  }

  // Check both arguments have length > 0
  if(seq_a[0] == '\0' || seq_b[0] == '\0')
  {
    fprintf(stderr, "Error: Sequences must have length > 0\n");
    fflush(stderr);

    if(cmd->print_fasta && seq_a_name != NULL && seq_b_name != NULL)
    {
      fprintf(stderr, "%s\n%s\n", seq_a_name, seq_b_name);
    }

    fflush(stderr);

    return;
  }

  smith_waterman_align(seq_a, seq_b, &scoring, sw);

  aligner_t *aligner = smith_waterman_get_aligner(sw);
  size_t len_a = aligner->score_width-1, len_b = aligner->score_height-1;

  printf("== Alignment %zu lengths (%lu, %lu):\n", alignment_index, len_a, len_b);

  if(cmd->print_matrices)
  {
    alignment_print_matrices(aligner);
  }

  // seqA
  if(cmd->print_fasta && seq_a_name != NULL)
  {
    fputs(seq_a_name, stdout);
    putc('\n', stdout);
  }

  if(cmd->print_seq)
  {
    fputs(seq_a, stdout);
    putc('\n', stdout);
  }

  // seqB
  if(cmd->print_fasta && seq_b_name != NULL)
  {
    fputs(seq_b_name, stdout);
    putc('\n', stdout);
  }

  if(cmd->print_seq)
  {
    fputs(seq_b, stdout);
    putc('\n', stdout);
  }

  putc('\n', stdout);

  if(!cmd->min_score_set)
  {
    // If min_score hasn't been set, set a limit based on the lengths of seqs
    // or zero if we're running interactively
    cmd->min_score = wait_on_keystroke ? 0
                       : scoring.match * MAX2(0.2 * MIN2(len_a, len_b), 2);

    #ifdef SEQ_ALIGN_VERBOSE
    printf("min_score: %i\n", cmd->min_score);
    #endif
  }

  fflush(stdout);

  size_t hit_index = 0;


  printf("score: %i\n", aligner->max_score);

  fputs("==\n", stdout);
  fflush(stdout);

  // Increment sequence alignment counter
  alignment_index++;
}

void align_pair_from_file(read_t *read1, read_t *read2)
{
  align(read1->seq.b, read2->seq.b,
       (read1->name.end == 0 ? NULL : read1->name.b),
       (read2->name.end == 0 ? NULL : read2->name.b));
}

int main(int argc, char* argv[])
{
  #ifdef SEQ_ALIGN_VERBOSE
  printf("VERBOSE: on\n");
  #endif

  sw_set_default_scoring();
  cmd = cmdline_new(argc, argv, &scoring, SEQ_ALIGN_SW_CMD);

  // Align!
  sw = smith_waterman_new();

  if(cmd->seq1 != NULL)
  {
    // Align seq1 and seq2
    align(cmd->seq1, cmd->seq2, NULL, NULL);
  }

  // Align from files
  size_t i, num_of_file_pairs = cmdline_get_num_of_file_pairs(cmd);
  for(i = 0; i < num_of_file_pairs; i++)
  {
    const char *file1 = cmdline_get_file1(cmd, i);
    const char *file2 = cmdline_get_file2(cmd, i);
    if(file1 != NULL && *file1 == '\0' && file2 == NULL) {
      wait_on_keystroke = 1;
      file1 = "-";
    }
    align_from_file(file1, file2, &align_pair_from_file, !cmd->interactive);
  }

  // Free memory for storing alignment results
  smith_waterman_free(sw);

  cmdline_free(cmd);

  return EXIT_SUCCESS;
}
