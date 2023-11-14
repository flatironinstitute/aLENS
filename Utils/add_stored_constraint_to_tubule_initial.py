#!/usr/bin/env python

from pathlib import Path
from read_write_utils import read_sylinder_ascii_file, read_links_ascii_file, TriBendLink


syl_path = Path.cwd() / 'TubuleInitial_old.dat'
new_syl_path = Path.cwd() / 'TubuleInitial.dat'
syls = read_sylinder_ascii_file(syl_path)
links = read_links_ascii_file(syl_path)


with new_syl_path.open('w') as tf:
    tf.write('Initial tubule config file\n\n')
    for syl in syls:
        tf.write(syl.to_string())
    for l in links:
        tf.write(l.to_string())
