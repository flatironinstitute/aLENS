#!/usr/bin/env python

from pathlib import Path
from read_write_utils import read_links_ascii_file, TriBendLink


syl_path = Path('TubuleInitial.dat')
new_syl_path = Path('TubuleInitial_new.dat')
links = read_links_ascii_file(syl_path)

# create new tribend links
tribend_links = []
for link in links[:-2]:  # assumes only one filament
    tribend_links += [TriBendLink(
        f'T {link._next_id} {link._prev_id} {link._next_id+1}')]


with new_syl_path.open('w') as tf, syl_path.open('r') as sf:
    for line in sf:
        tf.write(line)
    for l in links:
        tf.write(l.to_string())
    for tbl in tribend_links:
        tf.write(tbl.to_string())
