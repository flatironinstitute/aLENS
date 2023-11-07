#!/usr/bin/env python

"""@package docstring
File: runlog_funcs.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import sys
import re
import time

import numpy as np

from datetime import datetime
from pathlib import Path


def get_walltime(log_path):
    """Uses the log file to calculate the total time the simulation took.
    This will not work for restarted simulations. Might want to fix that.

    @param log_path TODO
    @return: TODO

    """
    with open(log_path, 'r') as rlf:
        pattern = re.compile(r'\[(\d+-\d+-\d+\s\d+:\d+:\d+\.\d+)\]')
        line = rlf.readline()
        while not pattern.search(line):
            line = rlf.readline()
        start_wtime = pattern.search(line).group(0)

        for line in reversed(rlf.readlines()):
            if not pattern.search(line):
                continue
            end_wtime = pattern.search(line).group(0)
            break

    stripstr = '[%Y-%m-%d %H:%M:%S.%f]'
    end_dt = datetime.strptime(end_wtime, stripstr)
    start_dt = datetime.strptime(start_wtime, stripstr)
    return end_dt - start_dt


def get_wt_timestep(log_path):
    """Get the wall time for every time step.
    (Might want to add an option to coarse grain.)

    @param log_path Path to run.log
    @return: Times between timesteps

    """
    with open(log_path, 'r') as rlf:
        time_ptrn = re.compile(
            r'\[(\d+-\d+-\d+\s\d+:\d+:\d+\.\d+)\].*CurrentStep\W+(\d+)')
        stripstr = '%Y-%m-%d %H:%M:%S.%f'
        wt_arr = []
        for line in rlf:
            re_match = time_ptrn.search(line)
            if re_match:
                if not int(re_match.group(2)):
                    continue
                wt_arr += [
                    datetime.strptime(
                        time_ptrn.search(line).group(1),
                        stripstr)]

    fconv = np.vectorize(lambda x: x.total_seconds())
    wt_arr = np.array(wt_arr)
    return fconv(wt_arr[1:] - wt_arr[:-1])


def calc_timestep_stats(log_path):
    """TODO: Docstring for calc_timestep_states.

    @param log_path TODO
    @return: TODO

    """
    wt_step_arr = get_wt_timestep(log_path)
    tot_wt = wt_step_arr[-1] - wt_step_arr[0]
    return (float(wt_step_arr.mean()),
            float(np.median(wt_step_arr)),
            float(wt_step_arr.std()),
            float(np.max(wt_step_arr)))


##########################################
if __name__ == "__main__":
    total_time = get_walltime(Path(sys.argv[1]))
    print(total_time)
    print(f"Total time in second: {total_time.total_seconds()}")
