#!/usr/bin/env python

"""@package docstring
File: initialize_scan.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
from pathlib import Path
from shutil import copy, rmtree
import zipfile
import re
# import chi_pet
from chi_pet import ChiCreate
from chi_pet.Chi import parse_args, ChiMain


BASE_EQUIL_DIR = Path(
    '/home/alamson/DATA/Chromatin/22-08_aLc1_line200-6400_equil_runs')


def get_final_syl_file(run_dir):
    """Return the location of the final directory

    @param run_dir Run directory with a result directory in it
    @return: String of file that was read

    """
    # TODO Check to see if there is a zip file instead of a result directory
    #      This will determine if you want to use the regular path or not
    # Find the result#-# directory with the largest #
    result_dirs = [d for d in (run_dir / 'result').iterdir() if d.is_dir()]
    large_result_dir = max(
        result_dirs, key=lambda x: int(re.findall(r'[0-9]+', x.name)[-1]))
    # Within this folder, find the last SylinderAscii_#.dat
    large_syl_file = max(large_result_dir.glob('SylinderAscii_*.dat'),
                         key=lambda x: int(re.findall(r'[0-9]+', x.name)[0]))
    # Read file to a string and return
    return large_syl_file.read_text()


def main():
    """TODO: Docstring for main.
    @return: TODO

    """
    root_dir = Path.cwd()
    root_files = [f for f in root_dir.iterdir() if f.is_file()
                  and f.suffix != '.py']

    sim_dir = root_dir / 'simulations'
    # Always start out clean
    rmtree(sim_dir)
    sim_dir.mkdir()

    equil_dirs = [rdir for rdir in BASE_EQUIL_DIR.glob('line*')]
    for eq_dir in equil_dirs:
        # Make a new run simulation
        run_dir = sim_dir / eq_dir.name
        run_dir.mkdir()
        # Copy files from current directory into the new run dir
        for f in root_files:
            copy(f, run_dir)
        # Get initial condition files and put into new run directory
        init_str = get_final_syl_file(eq_dir)
        # Write new initial configuration files to run dir
        init_file = run_dir / 'TubuleInitial.dat'
        init_file.write_text(init_str)

        # Run Chi or Chi-based implementation on this new directory
        # Need to create an opts object to pass the creation script
        opts = parse_args()
        opts.workdir = run_dir
        opts.create = [f for f in run_dir.iterdir() if f.is_file()
                       and f.suffix == '.yaml']
        opts.args_file = 'args.yaml'
        opts.non_yaml = [run_dir / 'aLENS.X', init_file]

        c = ChiMain(opts)


##########################################
if __name__ == "__main__":
    main()
