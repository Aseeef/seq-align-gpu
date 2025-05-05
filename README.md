# Smith-Waterman High Performance CPU Implementation

This project is a high-performance, parallel and vectorized implementation of the **Smith-Waterman algorithm** with **affine gap penalties**, designed specifically for **protein sequence alignment** using **AVX2 vectorization** and **OpenMP** multithreading.

## Features

* Full Smith-Waterman with affine gap scoring
* Configurable substitution matrix
* One-to-many alignment (query vs. database)
* Linear-space memory implementation
* AVX2 SIMD vectorization (16-bit and 32-bit paths)
* Multi-threaded with OpenMP
* Memory and cache optimizations

## Usage

### Prerequisites
* Implementation was developed and tested against x86 architecture on Linux. Support of other operating systems and architectures is not guaranteed.
* Must have build tools installed.
* Provided query and database files must be in FASTA format.

```bash
make            # builds the project
# the first parameter in --files is always the query sequence file and the second is always the database you are querying
bin/smith_waterman --substitution_matrix scoring/PAM250.txt --printfasta --files database/query.fasta database/database.fasta
```

## Repository Structure

* `src/` - main source code
* `test/` - tests for correctness
* `benchmarks/` - benchmarking utilities
* `Final Report.pdf` - project report and benchmark figures

## Notes

* Designed for protein sequences only
* Assumes amino acid alphabet maps to 32-character substitution matrix
* See Final Report.pdf for full technical details and benchmark results

## Author

This project was developed by **Muhammad Aseef Imran**, for high-performance bioinformatics alignment.


### Other Credits

This project builds upon and was originally based on Isaac Turner's implementation as found at:
[https://github.com/noporpoise/seq-align](https://github.com/noporpoise/seq-align)
