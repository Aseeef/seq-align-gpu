/*
 tests/tests.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Feb 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <time.h> // needed for rand()
#include <unistd.h>  // need for getpid() for getting setting rand number


#include "smith_waterman.h"

//
// Tests
//
const char *suite_name = NULL;
char suite_pass = 0;
int suites_run = 0, suites_failed = 0, suites_empty = 0;
int tests_in_suite = 0, tests_run = 0, tests_failed = 0;

#define QUOTE(str) #str
#define ASSERT(x) {tests_run++; tests_in_suite++; if(!(x)) \
  { fprintf(stderr, "failed assert [%s:%i] %s\n", __FILE__, __LINE__, QUOTE(x)); \
    suite_pass = 0; tests_failed++; }}

void SUITE_START(const char *name)
{
  suite_pass = 1;
  suite_name = name;
  suites_run++;
  tests_in_suite = 0;
}

void SUITE_END()
{
  printf("Testing %s ", suite_name);
  size_t suite_i;
  for(suite_i = strlen(suite_name); suite_i < 80-8-5; suite_i++) printf(".");
  printf("%s\n", suite_pass ? " pass" : " fail");
  if(!suite_pass) suites_failed++;
  if(!tests_in_suite) suites_empty++;
}

//
// Useful MACROs

#define die(fmt,...) do { \
  fprintf(stderr, "[%s:%i] Error: %s() "fmt"\n", __FILE__, __LINE__, __func__, ##__VA_ARGS__); \
  exit(EXIT_FAILURE); \
} while(0);

#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MIN(x,y) ((x) <= (y) ? (x) : (y))

// size is buffer size including NUL byte
static char* make_rand_seq(char *seq, size_t size)
{
  if(size == 0) { return seq; }
  size_t i, len = rand() % (size-1);
  const char bases[] = "acgt";
  for(i = 0; i < len; i++) seq[i] = bases[rand() & 3];
  seq[len] = '\0';
  return seq;
}

void sw_test_no_gaps_smith_waterman()
{
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *result = alignment_create(256);

  const char* seq_a = "gacag";
  const char* seq_b = "tgaagt";

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;

  bool case_sensitive = true;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend, case_sensitive);

  smith_waterman_align(seq_a, seq_b, &scoring, sw);

  smith_waterman_fetch(sw, result);
  ASSERT(strcmp(result->result_a, "ga") == 0 &&
         strcmp(result->result_b, "ga") == 0);

  smith_waterman_fetch(sw, result);
  ASSERT(strcmp(result->result_a, "ag") == 0 &&
         strcmp(result->result_b, "ag") == 0);

  alignment_free(result);
  smith_waterman_free(sw);
}


void test_sw()
{
  SUITE_START("Smith-Waterman");

  sw_test_no_gaps_smith_waterman();

  SUITE_END();
}

int main(int argc, char **argv)
{
  if(argc != 1)
  {
    printf("  Unused args '%s..'\n", argv[1]);
    printf("Usage: ./bit_array_test\n");
    exit(EXIT_FAILURE);
  }

  // Initialise random number generator
  srand((unsigned int)time(NULL) + getpid());

  printf("  Test seq-align C library:\n\n");

  // Test suites go here
  test_sw();

  printf("\n");
  printf(" %i / %i suites failed\n", suites_failed, suites_run);
  printf(" %i / %i suites empty\n", suites_empty, suites_run);
  printf(" %i / %i tests failed\n", tests_failed, tests_run);

  printf("\n THE END.\n");

  return tests_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
