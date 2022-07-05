"""@package docstring
File: gen_hydro_test.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""
# %%
from pathlib import Path
from time import time
import numpy as np
import yaml
from read_write_utils import read_sylinder_ascii_file, read_links_ascii_file
from hydro_utils import rotne_prager_tensor
# %%


def get_pos_vec_from_particles(particles):

    rvec = np.zeros((len(particles), 3))

    for i, p in enumerate(particles):
        rvec[i, :] = p.get_com()

    return rvec


def force_vec_from_links(pos_vec, links, spring_const, rest_length):
    r_vectors = pos_vec.reshape((pos_vec.size // 3, 3))
    force_vec = np.zeros(r_vectors.shape)
    for link in links:
        dr_vec = r_vectors[link._prev_id] - r_vectors[link._next_id]
        dr = np.linalg.norm(dr_vec)
        if dr > 0.:
            force_vec[link._prev_id] -= spring_const * (
                1. - (rest_length/dr)) * dr_vec
            force_vec[link._next_id] -= force_vec[link._prev_id]
    return force_vec.flatten()


run_path = Path('.')
# Read in run_config_file
with open(run_path / 'RunConfig.yaml', 'r') as yf:
    run_params = yaml.safe_load(yf)

# Read in initial data
particles = read_sylinder_ascii_file(run_path / "TubuleInitial.dat")
links = read_links_ascii_file(run_path / "TubuleInitial.dat")

# Put position vectors into array
pos_vec = get_pos_vec_from_particles(particles).flatten()
# Initialize forces
radius = particles[0].radius
visc = run_params['viscosity']
# Calculate mobility coefficient
mob_tensor = rotne_prager_tensor(pos_vec, visc, radius)
# Get spring parameters
spring_const = run_params['linkKappa']
rest_length = run_params['linkGap'] + 2. * radius
# Evolve system according to parameters
dt = run_params['dt']
steps = int(run_params['timeTotal']/dt)
for i in range(steps):
    # Get force vector from springs
    force_vec = force_vec_from_links(pos_vec, links, spring_const, rest_length)

    # Add force
    pos_vec += np.einsum('ij,j->i', mob_tensor, force_vec)*dt
    mob_tensor = rotne_prager_tensor(pos_vec, visc, radius)
pos_vec = pos_vec.reshape((pos_vec.size // 3, 3))

# Output data files
np.savetxt(run_path / 'TubuleFinal.dat', pos_vec)
