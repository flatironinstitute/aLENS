

"""@package docstring
File: gen_long_filament_from_data.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

import yaml
import argparse
import random
from pathlib import Path
import shutil
import numpy as np
from copy import deepcopy
from read_write_utils import (Sylinder, Protein, ProteinAscii,
                              get_sorted_dat_paths,
                              read_sylinder_ascii_file,
                              gen_id)
from gen_long_filament_from_data import gen_long_filament_from_data, parse_args
from gen_flexible_init import Links
from typing import Sequence, Any, Optional


def main(opts):
    nbeads = (200*np.power(2, np.arange(6))).tolist()
    paths = [Path.cwd() / f'line{i}' for i in nbeads]
    for path in paths:
        path.mkdir(exist_ok=True)

    files = ['RunConfig.yaml',
             'ProteinConfig.yaml',
             'aLENS.X',
             'gitversion.txt']
    for f in files:
        for d in paths:
            shutil.copy(Path.cwd() / f, d / f)

    for n, p in zip(nbeads, paths):
        opts.nbeads = n
        opts.output_dir = p
        gen_long_filament_from_data(opts)


##########################################
if __name__ == "__main__":
    opts = parse_args()
    main(opts)
