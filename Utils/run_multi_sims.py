#!/usr/bin/env python

"""@package docstring
File: run_multi_sims.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""


import os
from pathlib import Path
import subprocess

PARTS_PER_THREAD = 2000


def get_part_num(init_file):
    """TODO: Docstring for get_part_num.

    @param init_file TODO
    @return: TODO

    """
    count = 0
    with open(init_file, 'r') as f:
        for line in f:
            if line[0] == 'C' or line[0] == 'S':
                count += 1
    return count


def main():
    """TODO: Docstring for main.
    @return: TODO

    """
    root_dir = Path.cwd().resolve()
    # dir_paths = [d for d in next(Path.cwd().glob('line*'))]
    for d in root_dir.glob('line*'):
        os.chdir(d)
        init_file = d / 'TubuleInitial.dat'
        if not init_file.exists():
            raise IOError(f'File not found in {str(d)}.')

        part_num = get_part_num(init_file)
        thread_num = int(part_num / PARTS_PER_THREAD)
        thread_num = 1 if thread_num == 0 else thread_num
        os.environ['OMP_NUM_THREADS'] = str(thread_num)
        outfile = d / 'runlog.out'
        errfile = d / 'errlog.err'
        olog = outfile.open('w')
        elog = errfile.open('w')
        p = subprocess.Popen(['aLENS.X'], stdout=olog, stderr=elog, shell=True)

        os.chdir(root_dir)


##########################################
if __name__ == "__main__":
    main()
