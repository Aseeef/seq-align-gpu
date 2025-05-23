/*
 alignment_scoring_load.c
 url: https://github.com/noporpoise/seq-align
 maintainer: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain, no warranty
 date: Nov 2013
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h> // tolower isspace

#include "string_buffer/string_buffer.h"

#include "alignment_cmdline.h"
#include "alignment_scoring_load.h"

/*
 * Helper function to handle errors during matrix or substitution pair loading.
 * Exits the program after printing the error message, file path, and line number.
 *
 * Parameters:
 *   err_msg     - The error message to display.
 *   file_path   - The path of the file being processed.
 *   line_num    - The line number in the file where the error occurred.
 *   is_matrix   - Flag indicating whether the error occurred in a substitution matrix (1) 
 *                 or a substitution pair (0).
 */
static void _loading_error(const char *err_msg, const char *file_path,
                           int line_num, char is_matrix)
__attribute__((noreturn));

static void _loading_error(const char *err_msg, const char *file_path,
                           int line_num, char is_matrix) {
    if (is_matrix) fprintf(stderr, "Error: substitution matrix : %s\n", err_msg);
    else fprintf(stderr, "Error: substitution pairs : %s\n", err_msg);

    if (file_path != NULL) fprintf(stderr, "File: %s\n", file_path);
    if (line_num != -1) fprintf(stderr, "Line: %s\n", file_path);

    exit(EXIT_FAILURE);
}

/*
 * Function to load a scoring matrix for sequence alignment.
 * Reads a file and populates scoring values for character pair mutations based on
 * either a whitespace-separated or custom-separator format.
 *
 * Parameters:
 *   file           - A gzFile pointer to the input file containing the scoring matrix.
 *   file_path      - The file path of the input file (for error reporting).
 *   scoring        - Pointer to the scoring_t structure where loaded scores will be stored.
 *   case_sensitive - Flag indicating if characters should be treated case-sensitively 
 *                    (non-zero) or case-insensitively (0).
 */
void align_scoring_load_matrix(gzFile file, const char *file_path,
                               scoring_t *scoring, char case_sensitive) {
    StrBuf *sbuf = strbuf_new(500);
    size_t read_length;
    int line_num = 0;

    // Read first line (column headings)
    while ((read_length = strbuf_reset_gzreadline(sbuf, file)) > 0) {
        strbuf_chomp(sbuf);

        if (sbuf->end > 0 && sbuf->b[0] != '#' && // line is not empty, not comment
            !string_is_all_whitespace(sbuf->b)) // and not whitespace
        {
            // Read first line

            if (sbuf->end < 2) {
                _loading_error("Too few column headings", file_path, line_num, 1);
            }

            break;
        }

        line_num++;
    }

    if (line_num == 0 && sbuf->end <= 0) {
        _loading_error("Empty file", file_path, -1, 0);
    }

    // If the separator character is whitespace,
    // the set of whitespace characters is used
    char sep = sbuf->b[0];

    if ((sep >= (int) '0' && sep <= (int) '9') || sep == '-') {
        _loading_error("Numbers (0-9) and dashes (-) do not make good separators",
                       file_path, line_num, 0);
    }

    char *characters = (char *) malloc(sbuf->end);
    int num_of_chars = 0;

    if (isspace(sep)) {
        char *next = sbuf->b;

        while ((next = string_next_nonwhitespace(next + 1)) != NULL) {
            characters[num_of_chars++] = case_sensitive ? *next : tolower(*next);
        }

        // Now read lines below
        while ((read_length = strbuf_reset_gzreadline(sbuf, file)) > 0) {
            strbuf_chomp(sbuf);

            char *from_char_pos = string_next_nonwhitespace(sbuf->b);

            if (from_char_pos == NULL || sbuf->b[0] == '#') {
                // skip this line
                continue;
            }

            char from_char = case_sensitive ? *from_char_pos : tolower(*from_char_pos);
            char to_char;

            char *score_txt = sbuf->b + 1;
            int score;

            int i;
            for (i = 0; i < num_of_chars; i++) {
                to_char = characters[i];

                if (!isspace(*score_txt)) {
                    _loading_error("Expected whitespace between elements - found character",
                                   file_path, line_num, 1);
                }

                score_txt = string_next_nonwhitespace(score_txt + 1);

                char *strtol_last_char_ptr = score_txt;
                score = (int) strtol(strtol_last_char_ptr, &strtol_last_char_ptr, 10);

                // If pointer to end of number string hasn't moved -> error
                if (strtol_last_char_ptr == score_txt) {
                    _loading_error("Missing number value on line", file_path, line_num, 1);
                }

                scoring_add_mutation(scoring, from_char, to_char, score);

                score_txt = strtol_last_char_ptr;
            }

            if (*score_txt != '\0' && !string_is_all_whitespace(score_txt)) {
                _loading_error("Too many columns on row", file_path, line_num, 1);
            }

            line_num++;
        }
    } else {
        size_t i;

        for (i = 0; i < sbuf->end; i += 2) {
            if (sbuf->b[i] != sep) {
                _loading_error("Separator missing from line", file_path, line_num, 1);
            }

            char c = case_sensitive ? sbuf->b[i + 1] : tolower(sbuf->b[i + 1]);
            characters[num_of_chars++] = c;
        }

        int score;

        // Read rows
        while ((read_length = strbuf_reset_gzreadline(sbuf, file)) > 0) {
            strbuf_chomp(sbuf);

            char from_char = case_sensitive ? sbuf->b[0] : tolower(sbuf->b[0]);

            if (from_char == '#' || string_is_all_whitespace(sbuf->b)) {
                // skip this line
                continue;
            }

            char *str_pos = sbuf->b;

            int to_char_index = 0;
            char to_char;

            while (*str_pos != '\0') {
                to_char = characters[to_char_index++];

                if (*str_pos != sep) {
                    _loading_error("Separator missing from line", file_path, line_num, 1);
                }

                // Move past separator
                str_pos++;

                char *after_num_str = str_pos;
                score = (int) strtol(str_pos, &after_num_str, 10);

                // If pointer to end of number string hasn't moved -> error
                if (str_pos == after_num_str) {
                    _loading_error("Missing number value on line", file_path, line_num, 1);
                }

                if (to_char_index >= num_of_chars) {
                    _loading_error("Too many columns on row", file_path, line_num, 1);
                }

                scoring_add_mutation(scoring, from_char, to_char, score);

                str_pos = after_num_str;
            }

            line_num++;
        }
    }

    free(characters);
    strbuf_free(sbuf);
}