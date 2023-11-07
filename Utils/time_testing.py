#!/usr/bin/env python

"""@package docstring
File: time_testing.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import argparse
import os
from shutil import copy, rmtree
from pathlib import Path
from subprocess import run
import toml
import yaml
from pprint import pprint

from .runlog_funcs import get_walltime, get_wt_timestep, calc_timestep_stats


def parse_args():
    parser = argparse.ArgumentParser(
        prog='aa_controller.py',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-p", "--path", default=".",
                        help="Path used in aLENS Analysis functions.")

    parser.add_argument("-nt", "--omp_num_threads",
                        type=int,
                        default=0,
                        help="Run aLENS in a directory one lower, collect the runtime statistics and put them in file located in the analysis directory. Takes in number of steps to run.")

    parser.add_argument("-f ", "--force", action='store_true',
                        help="Force analysis to occur. Overwrite previous analysis done.")

    parser.add_argument("-v ", "--verbose", action='store_true',
                        help="Output more information to stdout.")

    parser.add_argument("-c ", "--clean", action='store_true',
                        help="Remove time_testing directory after finishing time analysis.")

    opts = parser.parse_args()

    # Post parsing changes to options
    opts.path = Path(opts.path).resolve()
    print(opts.path)

    opts.result_dir = opts.path / 'result'
    opts.analysis_dir = opts.path / 'analysis'
    opts.analysis_dir.mkdir(exist_ok=True)

    return opts


def run_time_testing(n_time_steps, opts):
    """TODO: Docstring for run_time_testing.
    @return: TODO

    """
    start_dir = Path.cwd()
    # Create mock directory
    tt_path = opts.path / 'time_test'
    if tt_path.exists():
        # Start with a clean directory
        rmtree(tt_path)
    tt_path.mkdir()
    files = [p for p in opts.path.glob('*')
             if p.is_file() and p.suffix != '.zip']
    for f in files:
        copy(f, tt_path)

    try:
        # cd into mock directory to run simulations
        os.chdir(tt_path)

        # Change parameter files
        rc_path = tt_path / 'RunConfig.yaml'
        with rc_path.open('r') as rcf:
            params = yaml.safe_load(rcf)
            omp_num_threads = params.get('omp_num_threads', 1)
            if opts.omp_num_threads:
                print("   ! Overwriting OMP_NUM_THREADS in RunConfig.yaml. !")
                omp_num_threads = opts.omp_num_threads
            print("Number of OMP threads =", omp_num_threads)
            dt = float(params['dt'])
            params['timeTotal'] = dt * float(n_time_steps)
            # Don't write out during this process
            params['timeSnap'] = 10 * dt * float(n_time_steps)
        with rc_path.open('w') as rcf:
            yaml.dump(params, rcf)

        # Run alens (with specified number of
        out_path = (tt_path / 'runlog.out')
        err_path = (tt_path / 'runlog.err')
        run(['./aLENS.X'],
            stdout=out_path.open('w'),
            stderr=err_path.open('w'),
            env=dict(OMP_NUM_THREADS=str(omp_num_threads), **os.environ),
            )

        # Analyze runlog.out for run information
        tot_walltime = get_walltime(out_path)

        stats = {'total walltime': float(tot_walltime.total_seconds())}
        stats['mean step walltime'], stats['median step walltime'], stats['std of step walltime'], stats['max step walltime'] = calc_timestep_stats(
            out_path)
        # Put time analysis file in the analysis directory
        # (create if necessary)
        analysis_dir = start_dir / 'analysis'
        analysis_dir.mkdir(exist_ok=True)
        with (analysis_dir / f'timing_OMP{omp_num_threads}_Nsteps{n_time_steps}.toml').open('w') as tf:
            toml.dump(stats, tf)
        if opts.verbose:
            print("Timing stats")
            pprint(stats, sort_dicts=False)

    except BaseException:
        raise
    finally:
        os.chdir(start_dir)
        if opts.clean:
            rmtree(tt_path)


##########################################
if __name__ == "__main__":
    print("Not implemented yet")
