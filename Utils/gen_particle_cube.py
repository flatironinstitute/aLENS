"""@package docstring
File: gen_particle_cube.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description: Generates 27 particles in a cube formation as an 
initial conditiion for testing 
"""

import numpy as np
from read_write_utils import Particle  # nopep8


CUBE_LENGTH = .5  # um
P_RADIUS = .01  # um
DIRECTOR = np.asarray([1., 0., 0.])
P_LENGTH = .000002  # um

particles = []
gid = 0
for x in range(3):
    for y in range(3):
        for z in range(3):
            if x == 1 and y == 1 and z == 1:
                continue
            start_pos = .5*CUBE_LENGTH * \
                (np.asarray([x, y, z]) - 1.)
            start_pos[0] -= P_LENGTH * .5
            particles += [Particle(start_pos, DIRECTOR,
                                   P_LENGTH, P_RADIUS, 0, gid)]
            gid += 1

particle_str = ""
with open('TubuleInitial.dat', 'w') as tf:
    tf.write("# Initial configuration of particles\n#\n")
    for p in particles:
        tf.write(p.to_string())
