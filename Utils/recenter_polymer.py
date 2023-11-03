#!/usr/bin/env python

"""@package docstring
File: gen_rigid_filament_protein_file.py
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

    parser = argparse.ArgumentParser(prog='recenter_polymer.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''
Generate a protein file that adds bending rigidity to simulations
Example parameter file with comments:
    TODO
'''
                                     )

    parser.add_argument("input",
                        help="SylinderAscii file to recenter an make a new TubuleInitial.dat.")
    parser.add_argument("-p", "--protein", default=None,
                        help="A ProteinAscii data file. This option will keep"
                        "all protein positions and states in this file.")
    opts = parser.parse_args()

    opts.input = Path(opts.input)
    if opts.protein:
        opts.protein = Path(opts.protein)
    if not opts.input.exists():
        raise IOError(
            " {} does not exist. Put in valid path.".format(opts.input))

    # with param_path.open('r') as pf:
    #     opts.params = yaml.safe_load(pf)

    return opts


def main(opts):
    proteins = []
    gen_gid = gen_id()

    # Read in all the data
    segments = read_sylinder_ascii_file(opts.input)
    links = read_links_ascii_file(opts.input)
    # Get the center of mass of all segments
    fil_com = np.array([seg.get_com() for seg in segments]).mean(axis=0)
    # Subtract COM from all segments to have system centered at (0,0,0)
    max_rad = 0.
    tot_seg_vol = 0.
    for seg in segments:
        seg.start_pos -= fil_com
        seg.end_pos -= fil_com
        # Find minimum radius sphere that contains all segments
        max_rad = max([np.linalg.norm(seg.start_pos),
                       np.linalg.norm(seg.end_pos),
                       max_rad])
        tot_seg_vol += (4./3.) * np.pi * (seg.radius**3)

    print(f"Min radius containing sphere: {max_rad}")
    bound_sphere_vol = (4./3.) * np.pi * max_rad**3
    print(f"Packing fraction: {tot_seg_vol/bound_sphere_vol}")

    with open("TubuleInitial.dat", 'w') as f:
        f.write('# Initial configuration of tubules\n#\n')
        for seg in segments:
            f.write(seg.to_string())
        for link in links:
            f.write(link.to_string())

    # Update proteins if protein file was given
    if opts.protein:
        proteins = read_protein_ascii_file(opts.protein)
        with open("ProteinInitial.dat", 'w') as f:
            f.write('# Initial configuration of proteins\n#\n')
            for prot in proteins:
                prot.end0_pos -= fil_com
                prot.end1_pos -= fil_com
                f.write(prot.to_string())


##########################################
if __name__ == "__main__":
    opts = parse_args()
    main(opts)
