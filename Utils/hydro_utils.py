
"""@package docstring
File: hydro_utils.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description: Python utilities for calculating hydrodynamic forces
"""
import numpy as np


def rpy_fluid_vel_mat(pos_vec, r_vec, f_vec, eta, a):
    '''
    The fluid velocity perturbation around particle from forces on fluid

    pos: Nx3 vector of particle positions
    r_vec: Mx3 vector of all the positions of all the vectors on the fluid
    f_vec: Mx3 vector of forces on fluid where N is the number of forces on fluid
    eta: viscosity of the fluid
    a: radius of particle i

    return vel: fluid velocity moving particle i
    '''

    factor = 1. / (8. * np.pi * eta)

    dr_vec = pos_vec[:, np.newaxis, :] - r_vec[np.newaxis, :, :]

    dr2 = np.einsum('ijk,ijk->ij', dr_vec, dr_vec)
    dr = np.sqrt(dr2)
    dr_dot_f = np.einsum('ijk,jk->ij', dr_vec, f_vec)

    # Off-diagnol term from mobility matrix
    oterm = np.divide(np.einsum('ij,ijk->ijk', dr_dot_f, dr_vec),
                      dr2[:, :, np.newaxis])

    term1 = np.divide(f_vec + oterm, dr[:, :, np.newaxis])
    term2 = np.divide(a*a*((f_vec/3.) - oterm), dr2[:, :, np.newaxis])
    return factor * np.einsum('ijk->ik', term1 + term2)


def rotne_prager_tensor(r_vectors, eta, a):
    '''
    Calculate free rotne prager tensor for particles at locations given by
    r_vectors of radius a.
    '''
    # Extract variables
    r_vectors = r_vectors.reshape((r_vectors.size // 3, 3))
    x = r_vectors[:, 0]
    y = r_vectors[:, 1]
    z = r_vectors[:, 2]
    # Compute distances between blobs
    dx = x - x[:, None]
    dy = y - y[:, None]
    dz = z - z[:, None]
    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    # Compute scalar functions f(r) and g(r)
    factor = 1.0 / (6.0 * np.pi * eta)
    fr = np.zeros_like(dr)
    gr = np.zeros_like(dr)
    sel = dr > 2.0 * a
    nsel = np.logical_not(sel)
    sel_zero = dr == 0.
    nsel[sel_zero] = False
    fr[sel] = factor * (0.75 / dr[sel] + a**2 / (2.0 * dr[sel]**3))
    gr[sel] = factor * (0.75 / dr[sel]**3 - 1.5 * a**2 / dr[sel]**5)
    fr[sel_zero] = (factor / a)
    fr[nsel] = factor * (1.0 / a - 0.28125 * dr[nsel] / a**2)
    gr[nsel] = factor * (3.0 / (32.0 * a**2 * dr[nsel]))
    # Build mobility matrix of size 3N \times 3N
    M = np.zeros((r_vectors.size, r_vectors.size))
    M[0::3, 0::3] = fr + gr * dx * dx
    M[0::3, 1::3] = gr * dx * dy
    M[0::3, 2::3] = gr * dx * dz
    M[1::3, 0::3] = gr * dy * dx
    M[1::3, 1::3] = fr + gr * dy * dy
    M[1::3, 2::3] = gr * dy * dz
    M[2::3, 0::3] = gr * dz * dx
    M[2::3, 1::3] = gr * dz * dy
    M[2::3, 2::3] = fr + gr * dz * dz

    return M
