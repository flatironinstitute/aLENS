#!/usr/bin/env python

"""@package docstring
File: make_lammps_from_alens.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

import yaml
import argparse
from pathlib import Path
from read_write_utils import (Sylinder, Protein, ProteinAscii,
                              read_sylinder_ascii_file,
                              read_protein_ascii_file,
                              read_links_ascii_file,
                              gen_id)
import numpy as np


def parse_args():

    parser = argparse.ArgumentParser(prog='make_lammps_from_alens.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''
Generate a lammps data file from an aLENS SylinderAscii or TubuleInitial.dat file
Example parameter file with comments:
    TODO
'''
                                     )

    parser.add_argument("input",
                        help=("SylinderAscii file or TubuleInitial.dat to"
                              " recenter an make a new TubuleInitial.dat."))
    parser.add_argument("-lj", "--use_lj_units", action="store_true",
                        help="Convert units to lennard-jones normalized units.")
    opts = parser.parse_args()

    opts.input = Path(opts.input)
    if not opts.input.exists():
        raise IOError(
            " {} does not exist. Put in valid path.".format(opts.input))

    return opts


def main(opts):
    # Read in all the data
    segments = read_sylinder_ascii_file(opts.input)
    links = read_links_ascii_file(opts.input)
    with (Path.cwd() / 'RunConfig.yaml').open('r') as rcf:
        run_params = yaml.safe_load(rcf)

    with open("data.filament", 'w') as f:
        f.write('# Initial LAMMPS filament configuration from alens files\n')

        # Write header for file
        f.write(f"""
{len(segments)} atoms
{len(links)} bonds
0 angles
0 dihedrals
0 impropers

1 atom types
1 bond types
0 angle types
0 dihedral types
0 improper types

""")

        scale_factor = 1. / (2. *
                             segments[0].radius) if opts.use_lj_units else 1.

        # Write dimensions of box
        sbhi = np.array(run_params['simBoxHigh']) * scale_factor
        sblo = np.array(run_params['simBoxLow']) * scale_factor

        f.write((f"{sblo[0]} {sbhi[0]} xlo xhi\n"
                 f"{sblo[1]} {sbhi[1]} ylo yhi\n"
                 f"{sblo[2]} {sbhi[2]} zlo zhi\n"))

        # Write extraneous header things
        f.write(f"""
Masses

1 {1.0 if opts.use_lj_units else .00000002}

""")

        # Write particle section
        f.write("Atoms\n\n")
        for seg in segments:
            f.write(seg.get_lammps_str_to_write(scale_factor=scale_factor))
        f.write("\n")

        # Write velocities
        # (start everything at zero since alens ascii files do not
        # store velocity data. We should equilibrate the system for
        # longer if possible.)
        f.write("Velocities\n\n")
        for seg in segments:
            f.write(f"{int(seg.gid) + 1} 0.0 0.0 0.0\n")
        f.write("\n")

        f.write("Bonds\n\n")
        bond_id_gen = gen_id(1)
        for link in links:
            f.write(link.get_lammps_str_to_write(next(bond_id_gen)))


##########################################
if __name__ == "__main__":
    opts = parse_args()
    main(opts)
