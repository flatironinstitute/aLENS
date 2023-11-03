#!/usr/bin/env python

"""@package docstring
File: gen_flexible_init.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import yaml
import argparse
from pathlib import Path
import numpy as np
from copy import deepcopy
from read_write_utils import gen_id, read_sylinder_ascii_file, Protein, Particle
from typing import Sequence, Any, Optional


GLOBAL_PARAM_KEYS = {
    "rng_seed": 1234,
    "data_file": "TubuleInitialOld.dat",
    "num_proteins": 0,
    "sys_length": 1.0,
}


def normalize(vl):
    v = np.array(vl)
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def set_default_parameters(params, default_params, check_for_none=False):
    """Function to set global parameters in

    @param params TODO
    @return: TODO

    """
    for key, val in default_params.items():
        params[key] = params.get(key, val)
        if check_for_none and params[key] is None:
            raise RuntimeError(
                f"Parameter '{key}' not initialized globally or locally.")


def parse_args():

    parser = argparse.ArgumentParser(prog='gen_spindle_config.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''

'''
                                     )
    parser.add_argument("input",
                        help="Initial condition parameter yaml file.")
    opts = parser.parse_args()

    param_path = Path(opts.input)
    if not param_path.exists():
        raise IOError(
            " {} does not exist. Put in valid path.".format(param_path))

    with param_path.open('r') as pf:
        opts.params = yaml.safe_load(pf)

    return opts


def main(opts):
    """Main loop for generating flexible filaments
    @return: TODO

    """
    ti_path = Path("TubuleInitialNew.dat")
    pi_path = Path("ProteinInitialNew.dat")
    # Initialize keys of global parameters
    set_default_parameters(opts.params, GLOBAL_PARAM_KEYS)

    with ti_path.open('w') as tif, pi_path.open('w') as pif:
        tif.write('# Initial configuration of sylinders and spheres \n#\n')

        # Check for old data file to copy over
        data_path = Path(opts.params['data_file'])
        sylinders = read_sylinder_ascii_file(data_path)
        max_syl_id = max([int(syl.gid) for syl in sylinders])
        syl_gen_id = gen_id(max_syl_id + 1)
        prot_gen_id = gen_id(0)

        # Write out sylinders to file
        for syl in sylinders:
            syl.grp_id = 0  # Re-write the group for all sylinders
            tif.write(syl.to_string())

        # Write out proteins (combination of spheres and crosslinkers)
        radius = 0.025
        length = 1e-6  # Just want this really small
        direction = [1, 0, 0]  # doesn't matter
        grp_id = 1
        for i in range(opts.params['num_proteins']):
            # Generate the body (sphere) of the protein randomly in the system
            # make a random numpy 3d vector
            pos = (np.random.rand(3)) * opts.params['sys_length']
            sphere = Particle(pos, direction, length, radius,
                              grp_id, next(syl_gen_id))
            tif.write(sphere.to_string())
            prot = Protein(next(prot_gen_id), 0,
                           sphere.start_pos, sphere.start_pos,
                           sphere.gid, -1)
            pif.write(prot.to_string())
            prot.gid = next(prot_gen_id)
            pif.write(prot.to_string())


##########################################
if __name__ == "__main__":
    opts = parse_args()
    main(opts)
