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
                              gen_id)


def parse_args():

    parser = argparse.ArgumentParser(prog='gen_rigid_filament_protein_file.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''
Generate a protein file that adds bending rigidity to simulations
Example parameter file with comments:
    TODO
'''
                                     )

    parser.add_argument("input",
                        help="Initial condition parameter yaml file.")
    parser.add_argument("-c", "--config", default=None,
                        help="A ProteinAscii data file. This option will keep"
                        "all protein positions and states in this file and "
                        "add new ones defined in input.")
    opts = parser.parse_args()

    param_path = Path(opts.input)
    if opts.config:
        opts.config = Path(opts.config)
    if not param_path.exists():
        raise IOError(
            " {} does not exist. Put in valid path.".format(param_path))

    with param_path.open('r') as pf:
        opts.params = yaml.safe_load(pf)

    return opts


def main(opts):
    fname = "ProteinInitial.dat"
    proteins = []
    gen_gid = gen_id()
    if opts.config:
        proteins = read_protein_ascii_file(opts.config)
        gen_gid = gen_id(proteins[-1].gid + 1)

    filaments = read_sylinder_ascii_file(Path('TubuleInitial.dat'))
    # Create proteins from position and GIDs of filaments
    for region in opts.params:
        fils = filaments[region['range'][0]:region['range'][1]]
        if region['step'] == 0:
            proteins += [Protein(next(gen_gid), region['type'],
                                 fil.start_pos[:], fil.start_pos[:],
                                 fil.gid, -1)
                         for fil in fils]
        else:
            proteins += [Protein(next(gen_gid), region['type'],
                                 fil0.start_pos[:], fil1.start_pos[:],
                                 fil0.gid, fil1.gid)
                         for fil0, fil1 in zip(fils[:-region['step']], fils[region['step']:])]
    with open(fname, 'w') as f:
        f.write('# Initial configuration of proteins\n#\n')
        for prot in proteins:
            f.write(prot.get_str_to_write())


##########################################
if __name__ == "__main__":
    opts = parse_args()
    main(opts)
