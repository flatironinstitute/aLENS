#!/usr/bin/env python3

import sys
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Utils.time_testing import run_time_testing   # nopep8

benchmark_root_dir = Path.cwd().resolve()

alens_exec_path = benchmark_root_dir / 'aLENS.X'
assert alens_exec_path.exists()

gitversion_path = benchmark_root_dir / 'gitversion.txt'
assert gitversion_path.exists()

# Make a dummy opts object to pass to run_time_testing


def opts(): return None


opts.clean = True
opts.verbose = False

omp_num_threads_list = [1, 2, 4, 8, 16]
# omp_num_threads_list = [4, 8]


# Aster gas ##########################################
os.chdir(benchmark_root_dir / 'AsterGas')
opts.path = Path.cwd()

ap = opts.path / 'aLENS.X'
if ap.exists():
    ap.unlink()
ap.symlink_to(alens_exec_path)

n_steps = 20
for omp_num_threads in omp_num_threads_list:
    opts.omp_num_threads = omp_num_threads
    run_time_testing(n_steps, opts)
    print(f"Aster gas {omp_num_threads}: Finished.")

print("Finished benchmarking aster gas.")


# Flexible condensed bead-spring filament ############
os.chdir(benchmark_root_dir / 'CollapsedStickyFilament')
opts.path = Path.cwd()

ap = opts.path / 'aLENS.X'
if ap.exists():
    ap.unlink()
ap.symlink_to(alens_exec_path)

n_steps = 20
for omp_num_threads in omp_num_threads_list:
    opts.omp_num_threads = omp_num_threads
    run_time_testing(n_steps, opts)

print("Finished benchmarking collapsed sticky filament.")

# Confined semi-flexible filaments ####################
os.chdir(benchmark_root_dir / 'ConfinedSemiFlexFilaments')
opts.path = Path.cwd()

ap = opts.path / 'aLENS.X'
if ap.exists():
    ap.unlink()
ap.symlink_to(alens_exec_path)

n_steps = 20
for omp_num_threads in omp_num_threads_list:
    opts.omp_num_threads = omp_num_threads
    run_time_testing(n_steps, opts)

print("Finished benchmarking confined semiflexible filaments.")
