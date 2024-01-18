#!/usr/bin/env python

"""@package docstring
File: gen_flexible_init.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""

from pathlib import Path
import numpy as np
from gen_flexible_init import FilSegment, Links


params = {
    'radius': .0035,
    'nsegs': 110,
    'seg_length': .035,
    'start_pos': np.array([0., 0., -1.93]),
    'director': np.array([0, 0, 1]),
    'grp_id': 0
}

segs = []
links = []
start_pos = params['start_pos']
for gid in range(params['nsegs']):
    segs += [FilSegment(start_pos,
                        params['director'],
                        params['seg_length'],
                        params['radius'],
                        params['grp_id'],
                        gid)]
    start_pos += params['seg_length'] * params['director']

    if gid > 0:
        links += [Links(gid-1, gid, 'P')]  # Pinned links
        links += [Links(gid-1, gid, 'B')]  # Bending constraint

ti_path = Path('TubuleInitial.dat')
pi_path = Path('ProtienInitial.dat')

with ti_path.open('w') as ti:
    ti.write('# Number of sylinders\n# {}\n'.format(len(segs)))
    for seg in segs:
        ti.write(seg.to_string())
    for link in links:
        ti.write(link.to_string())
