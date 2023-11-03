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
import numpy as np
from scipy.spatial.transform import Rotation
from read_write_utils import (Sylinder, Protein, ProteinAscii,
                              read_sylinder_ascii_file,
                              read_links_ascii_file,
                              gen_id)
from gen_flexible_init import Links


def normalize(vl):
    v = np.array(vl)
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


UNIT_CUBE_VERTEX_POSITIONS = np.array([[0., 0., 0.],
                                       [1., 0., 0.],
                                       [1., 1., 0.],
                                       [0., 1., 0.],
                                       [0., 1., 1.],
                                       [0., 0., 1.],
                                       [1., 0., 1.],
                                       [1., 1., 1.]])
ROT_UNIT_CUBE_VERT_POS = Rotation.from_rotvec(
    np.pi/4 * np.array([0., 0., 1.])).apply(UNIT_CUBE_VERTEX_POSITIONS)
ROT_CUBE_VERT_POS = Rotation.from_rotvec(
    np.arctan(np.sqrt(2)) * np.array([1., 0., 0.])).apply(ROT_UNIT_CUBE_VERT_POS)
# Normalize by the diagnol
ROT_CUBE_VERT_POS /= np.linalg.norm(ROT_CUBE_VERT_POS[-1])


def parse_args():

    parser = argparse.ArgumentParser(prog='gen_compact_fil_from_data.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''
Generate a sylinder file with 8x the number of beads that still maintains the same compacted shape and topology. 
'''
                                     )

    parser.add_argument("-i", "--input", default="TubuleInitial.dat",
                        help="A SylinderAscii data file that will be used to create a larger filament")

    opts = parser.parse_args()

    opts.syl_path = Path(opts.input)
    if not opts.syl_path.exists():
        raise IOError(
            " {} does not exist. Put in valid path.".format(opts.syl_path))
    return opts


def get_bead_cube_pos_from_orig_bead_info(pos, direction, radius):
    new_bead_radius = radius / (1. + np.sqrt(3.))
    vertex_pos = 2*(radius-new_bead_radius) * ROT_CUBE_VERT_POS

    # Find angles to rotate the cube
    azimuthal_angle = np.arctan2(direction[1], direction[0])
    polar_angle = np.arccos(direction[2])

    # Rotate cube
    # First add random rotation around the z-axis to break axial symmetry
    vertex_pos = Rotation.from_rotvec(
        np.random.rand() * 2. * np.pi * np.array([0., 0., 1.])).apply(vertex_pos)
    vertex_pos = Rotation.from_rotvec(
        polar_angle * np.array([0., 1., 0.])).apply(vertex_pos)
    vertex_pos = Rotation.from_rotvec(
        azimuthal_angle * np.array([0., 0., 1.])).apply(vertex_pos)

    # Translate cube
    vertex_pos += pos - (direction * (radius-new_bead_radius))

    return vertex_pos, new_bead_radius


def main(opts):
    gen_gid = gen_id()
    beads = read_sylinder_ascii_file(opts.syl_path)
    bead_radius = beads[0].radius
    # Spacing between beads
    seg_length = bead_radius * .0001

    # Todo: Find positions based off links. For now just assume a single filament
    # links = read_links_ascii_file(opts.syl_path)

    bead_pos_arr = []
    for ib in range(len(beads)-1):
        vec = normalize(beads[ib+1].start_pos - beads[ib].end_pos)
        pos = .5 * (beads[ib+1].start_pos + beads[ib].end_pos)

        cube_pos_arr, new_radius = get_bead_cube_pos_from_orig_bead_info(
            pos, vec, bead_radius)
        bead_pos_arr.append(cube_pos_arr)

    # Combine all bead positions into one array
    new_bead_pos_arr = np.concatenate(bead_pos_arr, axis=0)

    # Scale the positions of the beads
    new_bead_pos_arr *= bead_radius / new_radius

    # Generate new beads from combined positions
    new_beads = []
    for bead_pos in new_bead_pos_arr:
        start_pos = bead_pos - vec * seg_length
        end_pos = bead_pos + vec * seg_length
        new_beads.append(Sylinder("{0} {1} {2} {3:0.6f} {4:0.6f} {5:0.6f} {6:0.6f} {7:0.6f} {8:0.6f} {9}".format(
            "C",
            next(gen_gid), bead_radius,
            start_pos[0], start_pos[1], start_pos[2],
            end_pos[0], end_pos[1], end_pos[2],
            0)))

    # Add links between beads
    new_links = []
    for i in range(len(new_beads)-1):
        new_links.append(Links(new_beads[i].gid, new_beads[i+1].gid))

    # Write new beads and links to file
    with open(opts.syl_path.parent / "TubuleInitial_8x.dat", 'w') as f:
        for bead in new_beads:
            f.write(bead.to_string())
        for link in new_links:
            f.write(link.to_string())


##########################################
if __name__ == "__main__":
    opts = parse_args()
    main(opts)
