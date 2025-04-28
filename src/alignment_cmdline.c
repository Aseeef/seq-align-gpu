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
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <limits.h> // INT_MIN
#include <stdarg.h> // for va_list
#include <time.h>

#include "seq_file/seq_file.h"

#include "alignment.h"
#include "alignment_cmdline.h"

#include <smith_waterman.h>

#include "alignment_scoring_load.h"

char parse_entire_int(char *str, int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (int)tmp;
    return 1;
  }
}

char parse_entire_uint(char *str, unsigned int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(tmp > UINT_MAX || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (unsigned int)tmp;
    return 1;
  }
}

static void print_usage(enum SeqAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...)
  __attribute__((format(printf, 4, 5)))
  __attribute__((noreturn));

static void print_usage(enum SeqAlignCmdType cmd_type, score_t defaults[4],
                        const char *cmdstr, const char *errfmt, ...)
{
  if(errfmt != NULL) {
    fprintf(stderr, "Error: ");
    va_list argptr;
    va_start(argptr, errfmt);
    vfprintf(stderr, errfmt, argptr);
    va_end(argptr);

    if(errfmt[strlen(errfmt)-1] != '\n') fprintf(stderr, "\n");
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
"    --substitution_matrix <file>  see details for formatting\n\n",
          defaults[0], defaults[1],
          defaults[2], defaults[3]);

  if(cmd_type == SEQ_ALIGN_SW_CMD)
  {
    // SW specific
    fprintf(stderr,
"    --minscore <score>   Minimum required score\n"
"                         [default: match * MAX(0.2 * length, 2)]\n"
"    --maxhits <hits>     Maximum number of results per alignment\n"
"                         [default: no limit]\n"
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

void cmdline_free(cmdline_t *cmd)
{
  free(cmd);
}

#define usage(fmt,...) print_usage(cmd_type,defaults,argv[0],fmt, ##__VA_ARGS__)

cmdline_t* cmdline_new(int argc, char **argv, scoring_t *scoring,
                       enum SeqAlignCmdType cmd_type)
{
  cmdline_t* cmd = calloc(1, sizeof(cmdline_t));
  cmd->seq1 = cmd->seq2 = NULL;
  // All values initially 0

  // Store defaults
  score_t defaults[4] = {scoring->match, scoring->mismatch,
                         scoring->gap_open, scoring->gap_extend};

  if(argc == 1) usage(NULL);

  // First run through arguments to set up case_sensitive and scoring system

  // case sensitive needs to be dealt with first
  // (is is used to construct hash table for swap_scores)
  char substitutions_set = 0, match_set = 0, mismatch_set = 0;

  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(strcasecmp(argv[argi], "--help") == 0 ||
       strcasecmp(argv[argi], "-help") == 0 ||
       strcasecmp(argv[argi], "-h") == 0)
    {
      usage(NULL);
    }
    else if(strcasecmp(argv[argi], "--case_sensitive") == 0)
    {
      cmd->case_sensitive = 1;
    }
  }

  for(argi = 1; argi < argc; argi++)
  {
    if(argv[argi][0] == '-')
    {
        if(strcasecmp(argv[argi], "--case_sensitive") == 0)
      {
        // Already dealt with
        //case_sensitive = true;
      }
      else if(strcasecmp(argv[argi], "--printseq") == 0)
      {
        if(cmd_type != SEQ_ALIGN_SW_CMD)
          usage("--printseq only valid with Smith-Waterman");
        cmd->print_seq = true;
      }
      else if(strcasecmp(argv[argi], "--printmatrices") == 0)
      {
        cmd->print_matrices = true;
      }
      else if(strcasecmp(argv[argi], "--printfasta") == 0)
      {
        cmd->print_fasta = true;
      }
      else if(strcasecmp(argv[argi], "--pretty") == 0)
      {
        cmd->print_pretty = true;
      }
      else if(strcasecmp(argv[argi], "--colour") == 0)
      {
        cmd->print_colour = true;
      }
      else if(strcasecmp(argv[argi], "--stdin") == 0)
      {
        // Similar to --file argument below
        cmdline_set_files(cmd, "", NULL);
        cmd->interactive = true;
      }
      else if(argi == argc-1)
      {
        // All the remaining options take an extra argument
        usage("Unknown argument without parameter: %s", argv[argi]);
      }
      else if(strcasecmp(argv[argi], "--scoring") == 0)
      {
        // This handled above
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--substitution_matrix") == 0)
      {
        gzFile sub_matrix_file = gzopen(argv[argi+1], "r");
        if(sub_matrix_file == NULL) usage("Couldn't read: %s", argv[argi+1]);

        align_scoring_load_matrix(sub_matrix_file, argv[argi+1],
                                  scoring, cmd->case_sensitive);

        gzclose(sub_matrix_file);
        substitutions_set = true;

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--minscore") == 0)
      {
        if(cmd_type != SEQ_ALIGN_SW_CMD)
          usage("--minscore only valid with Smith-Waterman");

        if(!parse_entire_int(argv[argi+1], &cmd->min_score))
          usage("Invalid --minscore <score> argument (must be a +ve int)");

        cmd->min_score_set = true;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--maxhits") == 0)
      {
        if(cmd_type != SEQ_ALIGN_SW_CMD)
          usage("--maxhits only valid with Smith-Waterman");

        if(!parse_entire_uint(argv[argi+1], &cmd->max_hits_per_alignment))
          usage("Invalid --maxhits <hits> argument (must be a +ve int)");

        cmd->max_hits_per_alignment_set = true;

        argi++;
      }
      else if(strcasecmp(argv[argi], "--match") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->match))
        {
          usage("Invalid --match argument ('%s') must be an int", argv[argi+1]);
        }

        match_set = true;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--mismatch") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->mismatch))
        {
          usage("Invalid --mismatch argument ('%s') must be an int", argv[argi+1]);
        }

        mismatch_set = true;
        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapopen") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_open))
        {
          usage("Invalid --gapopen argument ('%s') must be an int", argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--gapextend") == 0)
      {
        if(!parse_entire_int(argv[argi+1], &scoring->gap_extend))
        {
          usage("Invalid --gapextend argument ('%s') must be an int",
                argv[argi+1]);
        }

        argi++; // took an argument
      }
      else if(strcasecmp(argv[argi], "--file") == 0)
      {
        cmdline_set_files(cmd, argv[argi+1], NULL);
        argi++; // took an argument
      }
      // Remaining options take two arguments but check themselves
      else if(strcasecmp(argv[argi], "--files") == 0)
      {
          printf("Adding %s and %s", argv[argi+1], argv[argi+2]);
        if(argi >= argc-2)
        {
          usage("--files option takes 2 arguments");
        }
        else if(strcmp(argv[argi+1], "-") == 0 && strcmp(argv[argi+2], "-") == 0)
        {
          // Read both from stdin
          cmdline_set_files(cmd, argv[argi+1], NULL);
        }
        else
        {
          cmdline_set_files(cmd, argv[argi+1], argv[argi+2]);
        }

        argi += 2; // took two arguments
      }
      else usage("Unknown argument '%s'", argv[argi]);
    }
    else
    {
      if(argc - argi != 2) usage("Unknown options: '%s'", argv[argi]);
      break;
    }
  }

  if(substitutions_set && !match_set)
  {
    // if substitution table set and not match/mismatch
    scoring->use_match_mismatch = 0;
  }

  if(scoring->use_match_mismatch && scoring->match < scoring->mismatch) {
    usage("Match value should not be less than mismatch penalty");
  }

  // Check for extra unused arguments
  // and set seq1 and seq2 if they have been passed
  if(argi < argc)
  {
    cmd->seq1 = argv[argi];
    cmd->seq2 = argv[argi+1];
  }

  if(cmd->seq1 == NULL && (cmd->file_path1 == NULL || cmd->file_path2 == NULL))
  {
    usage("No input specified");
  }

  return cmd;
}


void cmdline_set_files(cmdline_t *cmd, char* query, char* database)
{
  cmd->file_path1 = query;
  cmd->file_path2 = database;
}

char* cmdline_get_file1(cmdline_t *cmd)
{
  return cmd->file_path1;
}

char* cmdline_get_file2(cmdline_t *cmd)
{
  return cmd->file_path2;
}

static seq_file_t* open_seq_file(const char *path, bool use_zlib)
{
  return (strcmp(path,"-") != 0 || use_zlib) ? seq_open(path)
                                             : seq_dopen(fileno(stdin), false, false, 0);
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

void align_from_query_and_db(const char *query_path, const char *db_path, const scoring_t * scoring, sw_aligner_t * sw,
                     void (print_alignment)(const char *query_seq, const char *db_seq,
                                  const char *query_name, const char *db_name),
                     bool use_zlib)
{
    seq_file_t *query_file, *db_file;

    // Open query file
    if((query_file = open_seq_file(query_path, use_zlib)) == NULL)
    {
        fprintf(stderr, "Error: couldn't open query file %s\n", query_path);
        fflush(stderr);
        return;
    }

    // Open database file
    if((db_file = open_seq_file(db_path, use_zlib)) == NULL)
    {
        fprintf(stderr, "Error: couldn't open database file %s\n", db_path);
        fflush(stderr);
        seq_close(query_file);
        return;
    }

    // Read the single query sequence
    read_t query_read;
    seq_read_alloc(&query_read);

    if(seq_read(query_file, &query_read) <= 0)
    {
        fprintf(stderr, "Error: Query file %s is empty or invalid\n", query_path);
        fflush(stderr);
        seq_close(query_file);
        seq_close(db_file);
        seq_read_dealloc(&query_read);
        return;
    }

    struct timespec time_start, time_stop;
    double total_time = 0;
    // Read database sequences and align each with the query
    read_t db_read;
    seq_read_alloc(&db_read);

    while(seq_read(db_file, &db_read) > 0)
    {
        clock_gettime(CLOCK_REALTIME, &time_start);
        smith_waterman_align(query_read.seq.b, db_read.seq.b, scoring, sw);
        clock_gettime(CLOCK_REALTIME, &time_stop);
        total_time += interval(time_start, time_stop);
        print_alignment(query_read.seq.b, db_read.seq.b,
                (query_read.name.end == 0 ? NULL : query_read.name.b),
                (db_read.name.end == 0 ? NULL : db_read.name.b));
    }

    printf("Total time: %f\n", total_time);

    // Close files and free memory
    seq_close(query_file);
    seq_close(db_file);
    seq_read_dealloc(&query_read);
    seq_read_dealloc(&db_read);
}

// If seq2 is NULL, read pair of entries from first file
// Otherwise read an entry from each
void align_from_file(const char *path1, const char *path2,
                     void (align)(read_t *r1, read_t *r2),
                     bool use_zlib)
{
  seq_file_t *sf1, *sf2;

  if((sf1 = open_seq_file(path1, use_zlib)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path1);
    fflush(stderr);
    return;
  }

  if(path2 == NULL)
  {
    sf2 = sf1;
  }
  else if((sf2 = open_seq_file(path2, use_zlib)) == NULL)
  {
    fprintf(stderr, "Alignment Error: couldn't open file %s\n", path1);
    fflush(stderr);
    return;
  }

  // fprintf(stderr, "File buffer %zu zlib: %i\n", sf1->in.size, seq_use_gzip(sf1));

  read_t read1, read2;
  seq_read_alloc(&read1);
  seq_read_alloc(&read2);

  // Loop while we can read a sequence from the first file
  unsigned long alignments;

  for(alignments = 0; seq_read(sf1, &read1) > 0; alignments++)
  {
    if(seq_read(sf2, &read2) <= 0)
    {
      fprintf(stderr, "Alignment Error: Odd number of sequences - "
                      "I read in pairs!\n");
      fflush(stderr);
      break;
    }

    (align)(&read1, &read2);
  }

  // warn if no bases read
  if(alignments == 0)
  {
    fprintf(stderr, "Alignment Warning: empty input\n");
    fflush(stderr);
  }

  // Close files
  seq_close(sf1);

  if(path2 != NULL)
    seq_close(sf2);

  // Free memory
  seq_read_dealloc(&read1);
  seq_read_dealloc(&read2);
}
