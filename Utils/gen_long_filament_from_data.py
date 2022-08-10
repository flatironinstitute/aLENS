#!/usr/bin/env python

"""@package docstring
File: gen_long_filament_from_data.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import yaml
import argparse
import random
from pathlib import Path
import numpy as np
from copy import deepcopy
from read_write_utils import (Sylinder, Protein, ProteinAscii,
                              get_sorted_dat_paths,
                              read_sylinder_ascii_file,
                              gen_id)
from gen_flexible_init import Links
from typing import Sequence, Any, Optional


def parse_args():

    parser = argparse.ArgumentParser(prog='gen_long_filament_from_data.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''
    Generate an initialization file from equilibrated data file
    ''')
    parser.add_argument("input",
                        help="Result directory to take data files from.")
    parser.add_argument("-N", "--number", default=1, type=int,
                        help="Number of snapshots to use.")
    parser.add_argument("-n", "--nbeads", default=-1, type=int,
                        help="Number beads to use.")
    parser.add_argument("-s", "--start_xpos", default=0.,
                        help="Where on the x-axis to start.")
    parser.add_argument("-P", "--pinned", action='store_true',
                        help="Is this filament pinned at the ends.")
    parser.add_argument("-o", "--output_dir", default=Path.cwd(),
                        help="Where to put the InitialTubule.dat")
    opts = parser.parse_args()

    opts.input = Path(opts.input)
    if not opts.input.exists():
        raise IOError(
            " {} does not exist. Put in valid path.".format(opts.input.absolute()))
    return opts


def shift_filament(filaments, shift_vec):
    for fil in filaments:
        fil.start_pos += shift_vec
        fil.end_pos += shift_vec


def gen_long_filament_from_data(opts):
    gen_gid = gen_id()

    # Make a list of all the sylinder ascii files
    syl_files = get_sorted_dat_paths(opts.input)
    # Read the first file and get end to end distance
    filaments = read_sylinder_ascii_file(syl_files[0])
    fil_vec = filaments[-1].end_pos - filaments[0].start_pos
    shift_vec = np.array([opts.start_xpos, 0, 0]) - filaments[0].start_pos
    shift_direct = shift_vec/np.linalg.norm(shift_vec)

    # Loop over files chosen at random
    final_filament = []
    if opts.nbeads > 0:
        nbeads = opts.nbeads
        while nbeads > 0:
            filaments = read_sylinder_ascii_file(random.choice(syl_files))
            shift_filament(filaments, shift_vec)
            shift_vec += (fil_vec) + (filaments[0].radius*2.)*shift_direct
            for fil in filaments:
                final_filament += [fil]
                nbeads -= 1
                if nbeads <= 0:
                    break
    else:
        for i in range(opts.number):
            filaments = read_sylinder_ascii_file(random.choice(syl_files))
            shift_filament(filaments, shift_vec)
            shift_vec += (fil_vec) + (filaments[0].radius*2.)*shift_direct
            final_filament += filaments

    final_gid = len(final_filament)-1
    with open(opts.output_dir / "TubuleInitial.dat", 'w') as f:
        f.write('# Initial configuration of rods\n#\n')
        seg_str = ""
        link_str = ""
        for fil in final_filament:
            fil.gid = next(gen_gid)
            fil.type = 'S' if (
                (fil.gid == 0 or fil.gid == final_gid) and opts.pinned) else 'C'

            seg_str += fil.get_str_to_write()
        # Connect all the filament segments together
        links = [Links(i, i+1) for i in range(0, final_filament[-1].gid)]
        for link in links:
            link_str += link.get_str_to_write()
        f.write(seg_str)
        f.write(link_str)


##########################################
if __name__ == "__main__":
    opts = parse_args()
    gen_long_filament_from_data(opts)
