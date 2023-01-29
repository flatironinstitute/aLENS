#!/usr/bin/env python

"""@package docstring
File: make_last_step_new_init.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

from pathlib import Path


import yaml
import argparse
from shutil import copy, rmtree
import zipfile
import re
from glob import glob


def get_file_num(f):
    if isinstance(f, str):
        return int(re.findall(r'\d+', f)[-1])

    return int(re.findall(r'\d+', f.name)[-1])


def get_last_ascii_files(run_dir):
    """Return the location of the final directory

    @param run_dir Run directory with a result directory in it
    @return: String of file that was read

    """
    if (run_dir / 'result').exists():
        result_dir = run_dir / 'result'
        is_zip = False
    elif (run_dir / 'result.zip').exists():
        result_zip = zipfile.ZipFile(run_dir / 'result.zip')
        result_path = zipfile.Path(run_dir / 'result.zip')
        is_zip = True
    else:
        raise FileNotFoundError(
            f'Could not find result directory or zipfile in {run_dir}.')

    if is_zip:
        sy_reg = re.compile(r'.*SylinderAscii.*.dat')
        last_syl_file = result_path / max(list(filter(sy_reg.search, result_zip.namelist())),
                                          key=get_file_num)

        prot_reg = re.compile(r'.*ProteinAscii.*.dat')
        last_prot_file = result_path / max(list(filter(prot_reg.search, result_zip.namelist())),
                                           key=get_file_num)
    else:
        result_dirs = [d for d in result_dir.iterdir() if d.is_dir()
                       and d.name != 'PNG']
        large_result_dir = max(
            result_dirs, key=lambda x: int(re.findall(r'[0-9]+', x.name)[-1]))
        # Within this folder, find the last SylinderAscii_#.dat
        last_syl_file = max(large_result_dir.glob('SylinderAscii_*.dat'),
                            key=get_file_num)

        last_prot_file = max(large_result_dir.glob('ProteinAscii_*.dat'),
                             key=get_file_num)

    syl_snap_num = get_file_num(last_syl_file)
    prot_snap_num = get_file_num(last_prot_file)
    if syl_snap_num != prot_snap_num:
        raise AssertionError(
            'Snapshot numbers of sylinder and protein files do not match.')
    # Read file to a string and return
    return last_syl_file, last_prot_file


def modify_run_config(root_dir):
    """Modify the RunConfig file for time testing purposes.

    @param root_dir Directory with RunConfig.yaml in it
    @return: TODO

    """
    rc_path = root_dir / 'RunConfig.yaml'
    with rc_path.open('r') as rcy:
        params = yaml.safe_load(rcy)
        time_step = float(params['dt'])
        # If time testing, only run for 100 steps
        params['timeTotal'] = 100. * time_step
        # If time testing, make sure timeSnap is larger than totalTime
        # so no results are recorded, slowing down timing
        params['timeSnap'] = 10. * params['timeTotal']
        # Output maximum timing data to log files
        params['logLevel'] = 0
        params['timerLevel'] = 0
        params['initPreSteps'] = 0
    with rc_path.open('w') as rcf:
        yaml.dump(params, rcf)


def make_last_step_new_init(args):
    root_dir = args.root_dir
    last_syl_file, last_prot_file = get_last_ascii_files(root_dir)
    if last_syl_file == False:
        raise FileNotFoundError(
            'No result folder or file found. Cannot extract last files.')
        return

    tsi_file = Path('TimeStepInfo.txt')
    if tsi_file.exists() and not args.step_info:
        print('Removing TimeStepInfo.txt')
        tsi_file.unlink()

    if args.time_test:
        modify_run_config(root_dir)

    # Copy over files to initial configurations
    if isinstance(last_syl_file, zipfile.Path):
        with last_syl_file.open('r') as lsf, (root_dir / f'TubuleInitial{args.string_append}.dat').open('w') as tid:
            for line in lsf:
                tid.write(str(line, 'utf-8'))

        with last_prot_file.open('r') as lpf, (root_dir / f'ProteinInitial{args.string_append}.dat').open('w') as pid:
            for line in lpf:
                pid.write(str(line, 'utf-8'))
    else:
        copy(last_syl_file, root_dir / 'TubuleInitialNew.dat')
        copy(last_prot_file, root_dir / 'ProteinInitialNew.dat')


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-tt', '--time_test', action='store_true',
                        help='Modify RunConfig.yaml to only run for 100 steps without snapshots for quick timing testing.')
    parser.add_argument('-si', '--step_info', action='store_true',
                        help='Keep TimeStepInfo.txt file.')
    parser.add_argument('--string_append', type=str, default='',
                        help='Keep TimeStepInfo.txt file.')
    return parser.parse_args()


##########################################
if __name__ == "__main__":
    args = get_args()
    args.root_dir = Path.cwd()
    make_last_step_new_init(args)
