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
from read_write_utils import read_force_file, read_sylinder_ascii_file
from hydro_utils import rpy_fluid_vel_mat
# %%


def get_pos_vec_from_particles(particles):

    rvec = np.zeros((len(particles), 3))

    for i, p in enumerate(particles):
        rvec[i, :] = p.get_com()

    return rvec


def get_vecs_from_forces(forces):
    r_vec = np.zeros((len(forces), 3))
    f_vec = np.zeros((len(forces), 3))

    for i, f in enumerate(forces):
        r_vec[i, :] = f.pos[:]
        f_vec[i, :] = f.force[:]

    return r_vec, f_vec


run_path = Path('.')
# Read in run_config_file
with open(run_path / 'RunConfig.yaml', 'r') as yf:
    run_params = yaml.safe_load(yf)

# Read in initial data
particle_list = read_sylinder_ascii_file(run_path / "TubuleInitial.dat")
force_list = read_force_file(run_path / "FluidForce.dat")

# Put position vectors into array
pos_vec = get_pos_vec_from_particles(particle_list)
print(pos_vec)
# Initialize forces
r_vec, f_vec = get_vecs_from_forces(force_list)
visc = run_params['viscosity']
radius = particle_list[0].radius


vel_vec = rpy_fluid_vel_mat(pos_vec, r_vec, f_vec, visc, radius)

# Evolve system according to parameters
dt = run_params['dt']
steps = int(run_params['timeTotal']/dt)
for i in range(steps):
    pos_vec += vel_vec*dt
    vel_vec = rpy_fluid_vel_mat(pos_vec, r_vec, f_vec, visc, radius)

# Output data files
np.savetxt(run_path / 'TubuleFinal.dat', pos_vec)
