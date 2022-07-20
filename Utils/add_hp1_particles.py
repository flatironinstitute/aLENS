#!/usr/bin/env python

"""@package docstring
File: add_hp1_particles.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""


import yaml
import sys
import argparse
import random
from pathlib import Path
import numpy as np
from copy import deepcopy
from read_write_utils import (Sylinder, Protein, ProteinAscii,
                              get_sorted_dat_paths,
                              read_sylinder_ascii_file,
                              gen_id)
from gen_flexible_init import Links, FilSegment
from typing import Sequence, Any, Optional


#######################################################################
#                      Parameters to be changed                       #
#######################################################################

STICKY_RANGE = (7000, 9000)

CENTER_SHIFT = np.array([50., 0., 0.])
SEG_LENGTH = .00001
SEG_RADIUS = .01
HP1_NUM = 2000
STICKY_PER_HP1 = 2
RNG = np.random.default_rng()


def get_random_vec():
    """Get random unit vector
    @return: TODO

    """
    theta = np.arccos(RNG.uniform(-1, 1))
    phi = RNG.uniform(0, 2 * np.pi)
    return np.asarray([np.sin(theta) * np.cos(phi),
                       np.sin(theta) * np.sin(phi),
                       np.cos(theta)])


def make_sticky_particles(fil_id_gen, sys_rad):
    """Make spheres surrounding particles

    @param id_gen ID generator

    """
    prot_id_gen = gen_id()
    # sys_rad = float(self.kwargs['sys_dim'][0])
    sphrs = []
    prots = []
    for i in range(HP1_NUM):
        gid = next(fil_id_gen)
        director = get_random_vec()  # Put on surface of sphere to start
        end_pos = sys_rad * get_random_vec() + CENTER_SHIFT

        sphrs += [FilSegment(end_pos, director,
                             SEG_LENGTH, SEG_RADIUS,
                             1, gid)]
        for j in range(STICKY_PER_HP1):
            prots += [Protein(next(prot_id_gen), 0, end_pos, end_pos, gid, -1)]
    return sphrs, prots


def main():

    # Make a list of all the sylinder ascii files
    filaments = read_sylinder_ascii_file(Path(sys.argv[1]))

    seg_str = ""
    link_str = ""
    prot_str = ""
    for fil in filaments:
        # fil.type = 'S' if fil.gid == 0 or fil.gid == final_gid else 'C'
        fil.grp_id = 1 if int(
            fil.gid) >= STICKY_RANGE[0] and int(
            fil.gid) < STICKY_RANGE[1] else 0
        seg_str += fil.get_str_to_write()

    # Create crowding molecules
    gen_gid = gen_id(filaments[-1].gid + 1)
    sphrs, prots = make_sticky_particles(gen_gid, 10)
    for sphr in sphrs:
        seg_str += sphr.get_str_to_write()
    for prot in prots:
        prot_str += prot.get_str_to_write()
    # Connect all the filament segments together
    links = [Links(i, i + 1) for i in range(0, filaments[-1].gid)]
    for link in links:
        link_str += link.get_str_to_write()

    with open("TubuleInitial.dat", 'w') as f:
        f.write('# Initial configuration of rods\n#\n')
        f.write(seg_str)
        f.write(link_str)
    with open("ProteinInitial.dat", 'w') as f:
        f.write('# Initial configuration of proteins\n#\n')
        f.write(prot_str)


##########################################
if __name__ == "__main__":
    main()
