"""

This script was generated with AI Help.
Prompt: https://chatgpt.com/share/680ee10d-383c-800f-8374-07772828369c

"""

import subprocess
import re
import os
from tqdm import tqdm
from statistics import mean

# Configuration
executables_dir = '../executables'
command_args = ['--substitution_matrix', '../scoring/PAM250.txt', '--files', '../database/query.fasta', '../database/database.fasta']
ome_thresholds = [1, 2, 4, 8, 16, 32, 64]
repeats = 6

def run_benchmark(executable_path, omp_threads=None):
    times = []
    env = os.environ.copy()
    if omp_threads is not None:
        env['OMP_NUM_THREADS'] = str(omp_threads)

    for _ in tqdm(range(repeats), desc=executable_path.split("/")[-1] + ("" if omp_threads is None else f" (OMP_NUM_THREADS={omp_threads})")):
        try:
            result = subprocess.run([executable_path] + command_args, capture_output=True, text=True, env=env, check=True)
            match = re.search(r'Total time: ([0-9]*\.?[0-9]+)', result.stdout)
            if match:
                times.append(float(match.group(1)))
            else:
                print(f"Warning: 'Total time' not found in output of {executable_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error running {executable_path}: {e}")

    return mean(times) if times else None

def main():
    results = []
    for exe in os.listdir(executables_dir):
        exe_path = os.path.join(executables_dir, exe)

        if not os.path.isfile(exe_path) or not os.access(exe_path, os.X_OK):
            continue  # Skip non-executable files

        if 'omp' in exe:
            for threads in ome_thresholds:
                avg_time = run_benchmark(exe_path, omp_threads=threads)
                if avg_time is not None:
                    results.append((exe, threads, avg_time))
        else:
            avg_time = run_benchmark(exe_path)
            if avg_time is not None:
                results.append((exe, None, avg_time))

    print("\nBenchmark Results:")
    for exe, threads, avg_time in results:
        thread_info = f" (OMP_NUM_THREADS={threads})" if threads is not None else ""
        print(f"{exe}{thread_info}: Average Total Time = {avg_time:.6f} seconds")

if __name__ == "__main__":
    main()
