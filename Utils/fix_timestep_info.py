#!/usr/bin/env python

"""@package docstring
File: fix_timestep_info.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
from pathlib import Path
import yaml
from shutil import copy, rmtree
import zipfile
import re


def get_final_syl_file(run_dir):
    """Return the location of the final directory

    @param run_dir Run directory with a result directory in it
    @return: String of file that was read

    """
    result_dir = run_dir / 'result'
    if not result_dir.exists():
        return False

    # Find the result#-# directory with the largest #
    result_dirs = [d for d in result_dir.iterdir() if d.is_dir()]
    large_result_dir = max(
        result_dirs, key=lambda x: int(re.findall(r'[0-9]+', x.name)[-1]))
    # Within this folder, find the last SylinderAscii_#.dat
    large_syl_file = max(large_result_dir.glob('SylinderAscii_*.dat'),
                         key=lambda x: int(re.findall(r'[0-9]+', x.name)[0]))
    # Read file to a string and return
    return large_syl_file


def main():
    """TODO: Docstring for main.
    @return: TODO

    """
    root_dir = Path.cwd()
    final_file = get_final_syl_file(root_dir)
    tsi_file = Path('TimeStepInfo.txt')
    if final_file == False:
        print('No results file. Starting over.')
        if tsi_file.exists():
            tsi_file.unlink()
        return

    with Path('RunConfig.yaml').open('r') as rcy:
        params = yaml.safe_load(rcy)

    snap_num = int(re.findall(r'\d+', final_file.name)[0])
    pvtp_file = final_file.parent / f'Sylinder_{snap_num}.pvtp'
    assert pvtp_file.exists()
    with final_file.open('r') as ff:
        ff.readline()
        steps_taken = int(float(re.findall(r'\d*\.?\d+', ff.readline())[0]) /
                          float(params['dt']))
    tsi_file.write_text(f"""{params['rngSeed']}
{steps_taken}
{snap_num}
{pvtp_file.name}""")


##########################################
if __name__ == "__main__":
    main()
