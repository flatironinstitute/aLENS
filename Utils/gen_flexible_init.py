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
from read_write_utils import gen_id, read_sylinder_ascii_file, read_protein_ascii_file, read_links_ascii_file
from typing import Sequence, Any, Optional


GLOBAL_PARAM_KEYS = {
    'radius': None,
    'gen_type': None,
    'nsegs': None,
    'seg_length': None,
    'sys_dim': None,
    'start_pos': None,
    'director': None,
    'link_gap': None,
    'link_type': None
}


def normalize(vl):
    v = np.array(vl)
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def parse_args():

    parser = argparse.ArgumentParser(prog='gen_flexible_init.py',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description='''
Generate flexible filament configuation file from parameter yaml file.
Example parameter file with comments

global_parameters:      # All filaments will use these unless specified in filament list
    radius: .01         # 10 nm radius for filament
    sys_dim: [.5,.5,.5] # Dimensions of system
    link_gap: .0002     # Gap between sylinders connected by link (um)

filaments:
    - gen_type: line            # Type of generated lines options are:
                                #     line, hilbert, crowder, and random_walk
      nsegs: 100                # Number of segments in the line
      seg_length: .000001       # Length of segments.
      start_pos: [0., 0.01, 0.] # Starting position of line
      director: [1.,0.,0.]      # Director of the starting position
      stat_segs: [0]            # GIDs of segments that will not move
      # Connect this filament's last segment with the first of the next.
      join_next: True
      add_links: [[0,2]]        # List of tuples of extra linked segments
      group: 0                  # All segments will be given this group id
    - gen_type: line
      nsegs: 100
      seg_length: .000001
      start_pos: [1.0, -0.01, 0.]
      director: [-1.,0.,0.]
      add_links: [[50, 149], [49, 150], [48, 151], [47, 152]]
      join_next: False
      group: 0

'''
                                     )
    parser.add_argument("input",
                        help="Initial condition parameter yaml file.")
    parser.add_argument("-d", "--directory", type=Path, default=Path.cwd(),
                        help="Directory to add initial file to. "
                        "Default is current working directory")
    opts = parser.parse_args()

    param_path = Path(opts.input)
    if not param_path.exists():
        raise IOError(
            " {} does not exist. Put in valid path.".format(param_path))

    with param_path.open('r') as pf:
        opts.params = yaml.safe_load(pf)

    return opts


def hilbert_3d(s, i, pos_arr, cur_pos, dr1, dr2, dr3):
    """Recursive function to create a 3D hilbert curve given a number of lattice
    points per side

    @param s Number of lattice points per side
    @return: TODO

    """
    if s == 1:
        pos_arr[i, :] = cur_pos
        return i + 1

    s = np.floor(s / 2)
    assert (s >= 1)

    cur_pos_ = np.copy(cur_pos)
    dr1_ = np.copy(dr1)
    dr2_ = np.copy(dr2)
    dr3_ = np.copy(dr3)
    for dr in (dr1_, dr2_, dr3_):
        cur_pos_ -= s * (dr < 0).astype(int) * dr

    i = hilbert_3d(s, i, pos_arr, cur_pos_, dr2_, dr3_, dr1_)
    i = hilbert_3d(s, i, pos_arr, cur_pos_ + s * dr1_, dr3_, dr1_, dr2_)
    i = hilbert_3d(s, i, pos_arr, cur_pos_ + s * (dr1_ + dr2_),
                   dr3_, dr1_, dr2_)
    i = hilbert_3d(s, i, pos_arr, cur_pos_ + s * dr2_,
                   -dr1_, -dr2_, dr3_)
    i = hilbert_3d(s, i, pos_arr, cur_pos_ + s * (dr2_ + dr3_),
                   -dr1_, -dr2_, dr3_)
    i = hilbert_3d(s, i, pos_arr, cur_pos_ + s * (dr1_ + dr2_ + dr3_),
                   -dr3_, dr1_, -dr2_)
    i = hilbert_3d(s, i, pos_arr, cur_pos_ + s * (dr1_ + dr3_),
                   -dr3_, dr1_, -dr2_)
    i = hilbert_3d(s, i, pos_arr, cur_pos_ + s * dr3_,
                   dr2_, -dr3_, -dr1_)
    return i


def spiral_initial(start_pos: Sequence[float],
                   direction: Sequence[float],
                   seg_length: float,
                   n_segments: int,
                   end_sep: float = 5,  # Separation between start and end points
                   n_periods: int = 1,
                   **kwargs):
    """ Generate a flexible filament starting and ending at two points creating a spiral in between.
    """
    chain_length = seg_length * n_segments
    freq = np.pi * 2 * n_periods / (n_segments - 1)
    amp = np.sqrt((chain_length**2) - (end_sep**2)) / (2 * np.pi * n_periods)

    def get_pos(x): return start_pos + \
        np.array([float(end_sep * x / n_segments),
                  amp * np.cos(freq * x),
                  -amp * np.sin(freq * x)])

    pos_arr = []
    direct_arr = []
    for i in range(n_segments):
        pos_arr += [get_pos(i)]
        # direct_arr += [get_direct(i)]

    pos_arr = np.array(pos_arr)
    # Make direction array
    direct_arr = np.diff(pos_arr, axis=0)
    direct_arr /= np.linalg.norm(direct_arr, axis=1)[:, None]
    direct_arr = np.concatenate([direct_arr, [direct_arr[-1]]], axis=0)

    return pos_arr, direct_arr


class FilSegment():

    """Stores information about each filament segment including position, gid,
    filament group, and connecting segments. """

    def __init__(self, start_pos, director, length, radius, grp_id, gid,
                 seg_type='C'):
        self.director = normalize(deepcopy(director))
        self.length = length
        self.radius = radius
        self.start_pos = deepcopy(start_pos)
        self.end_pos = length * self.director + start_pos
        self.grp_id = grp_id
        self.gid = gid
        self.type = seg_type

    def to_string(self):
        return '{0} {1} {2} {3:0.6f} {4:0.6f} {5:0.6f} {6:0.6f} {7:0.6f} {8:0.6f} {9}\n'.format(
            self.type, self.gid, self.radius,
            self.start_pos[0], self.start_pos[1], self.start_pos[2],
            self.end_pos[0], self.end_pos[1], self.end_pos[2],
            self.grp_id)


class Links():

    """Spring-like links between filament"""

    def __init__(self, prev_id, next_id, link_type='E'):
        """Initialize with the ids of the two objects to connect

        @param prev_id First object to connect
        @param next_id Second object to connect

        """

        self._prev_id = prev_id
        self._next_id = next_id
        self._link_type = link_type

    def to_string(self):
        """Return string that defines link"""
        return f'{self._link_type} {self._prev_id} {self._next_id}\n'


class FlexFilament():
    def __init__(self, start_pos=[0., 0., 0.], director=[0., 0., 1.],
                 radius=.01, seg_length=.000001, nsegs=100,
                 group=0, gen_type='line', link_gap=.0002, **kwargs):
        self.start_pos = deepcopy(start_pos)
        self.end_pos = deepcopy(self.start_pos)
        self.director = normalize(director)
        self.seg_length = seg_length
        self.radius = radius
        self.epsilon = 2. * self.radius + link_gap
        self.nsegs = nsegs
        self.grp_id = group  # defined outisde the classv
        self.type = gen_type
        self.rng = np.random.default_rng()
        self.kwargs = kwargs
        self.stat_segs = kwargs.get('stat_segs', [])
        self.extra_links = kwargs.get('add_links', [])

        self.length = (self.seg_length + self.epsilon) * nsegs

        self.segs = []
        self.links = []

    def get_director_for_fil_type(self, i, director_arr=None):
        if self.type == "line":
            return self.director
        if self.type == "random_walk":
            return self.get_random_vec()
        if self.type == "hilbert" and i < self.nsegs - 1:
            return director_arr[i, :]
        if self.type == "spiral":
            return director_arr[i]
        else:
            raise RuntimeError("Unknown filament type")

    def make_flex_fil(self, id_gen):
        """Create a string of sylinder strings to write to file

        @param id_gen ID generator

        """
        prev_id = -1  # start without connection

        # Preallocate arrays
        director_arr = None
        if self.type == 'crowder':
            self.make_crowders(id_gen)
            return
        elif self.type == "hilbert":
            director_arr = self.get_hilbert_arr()
        elif self.type == "spiral":
            pos_arr, director_arr = spiral_initial(self.start_pos, self.director,
                                                   self.seg_length + self.epsilon,
                                                   self.nsegs, **self.kwargs)

        for i in range(self.nsegs):
            gid = next(id_gen)
            seg_type = 'S' if gid in self.stat_segs else 'C'

            # Cases for different types of filament generation schemes
            director = self.get_director_for_fil_type(i, director_arr)

            # Create new segment using the director and end position
            self.segs += [FilSegment(self.end_pos, director,
                                     self.seg_length, self.radius,
                                     self.grp_id, gid, seg_type)]

            if self.type == "spiral":
                self.end_pos = pos_arr[i]
            else:
                self.end_pos += (self.seg_length + self.epsilon) * director

            # Add links to list only if they are not the first or last object
            if prev_id >= 0 and i < self.nsegs:
                self.links += [Links(prev_id, gid)]

            prev_id = gid

        # Additional links
        for link in self.extra_links:
            self.links += [Links(link[0], link[1])]

        # Connect to the next filament. Will break simulations if there is no
        # next filament
        if self.kwargs.get('join_next', False):
            self.links += [Links(prev_id, prev_id + 1)]

        print(f'Start position: {self.start_pos}')
        self.end_pos = list(self.segs[-1].end_pos)
        print(
            f'End position: {self.end_pos}')

    def make_semiflexible_fil(self, id_gen):
        """Create a string of sylinder strings to write to file

        @param id_gen ID generator

        """
        prev_id = -1

        for i in range(self.nsegs):
            gid = next(id_gen)
            seg_type = 'S' if gid in self.stat_segs else 'C'

            # Cases for different types of filament generation schemes
            if self.type == "random_walk":
                director = self.get_random_vec()

            # Create new segment using the director and end position
            self.segs += [FilSegment(end_pos, director,
                                     self.seg_length, self.radius,
                                     self.grp_id, gid, seg_type)]

            # Add links to list only if they are not the first or last object
            if prev_id >= 0 and i < self.nsegs:
                self.links += [Links(prev_id, gid, self.kwargs['link_type'])]
            if prev_id >= 0 and i < self.nsegs:
                self.links += [Links(prev_id, gid, B)]

    def make_crowders(self, id_gen):
        """Make crowding particles inside of sphere

        @param id_gen ID generator

        """
        sys_rad = float(self.kwargs['sys_dim'][0])
        for i in range(self.nsegs):
            gid = next(id_gen)
            seg_type = 'S' if gid in self.stat_segs else 'C'
            director = self.get_random_vec()  # Put on surface of sphere to start
            start_pos = sys_rad * self.get_random_vec()

            self.segs += [FilSegment(start_pos, director,
                                     self.seg_length, self.radius,
                                     self.grp_id, gid, seg_type)]

    def get_fil_seg_str(self):
        """ Loop over segments and create initialization string
        @return: String of filament segments

        """
        seg_str = ""
        for seg in self.segs:
            seg_str += seg.to_string()
        return seg_str

    def get_fil_link_str(self):
        """ Loop over links and create initialization string
        @return: String of links

        """
        link_str = ""
        for link in self.links:
            link_str += link.to_string()
        return link_str

    def get_hilbert_arr(self):
        """Generate positions of segments in 3D space filling curve.

        @return: TODO

        """
        # Figure out the size of the cube based off the number of segments
        side_points = 2
        while self.nsegs > (side_points**3):
            side_points *= 2
        ind = 0
        # Make list of lists to hand to the recursive hilbert curve generator
        pos_arr = np.zeros((side_points**3, 3))
        # Create the first point and displacement vectors for program
        cur_pos = np.zeros(3)
        dr1, dr2, dr3 = [np.asarray((1, 0, 0)),
                         np.asarray((0, 1, 0)),
                         np.asarray((0, 0, 1))]
        # Generate list
        ind = hilbert_3d(side_points, ind, pos_arr, cur_pos, dr1, dr2, dr3)
        # Make director by looking at the next segment
        return pos_arr[1:, :] - pos_arr[:-1, :]

    def get_random_vec(self):
        """Get random unit vector
        @return: TODO

        """
        # TODO Make it so that random vector only samples hemisphere away from
        # previous segment
        theta = np.arccos(self.rng.uniform(-1, 1))
        phi = self.rng.uniform(0, 2 * np.pi)
        return np.asarray([np.sin(theta) * np.cos(phi),
                           np.sin(theta) * np.sin(phi),
                           np.cos(theta)])

    def get_next_fil_section_start(self, director):
        direct = normalize(director)
        return (self.segs[-1].end_pos
                + direct * (self.radius + self.epsilon))

    def __repr__(self):
        return "Filament()"

    def __str__(self):
        return 'Filament {0}:\n  pos_start: {1}\n  pos_end: {2}\n  length: {3}\n  radius: {4}'.format(
            self.gid, self.pos_start, self.end_pos, self.length, self.radius)


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


def gen_flexible_init(opts):
    """Main loop for generating flexible filaments
    @return: TODO

    """
    fname = opts.directory / "TubuleInitial.dat"
    # Initialize keys of global parameters
    set_default_parameters(opts.params['global_parameters'], GLOBAL_PARAM_KEYS)

    with fname.open('w') as f:
        f.write('# Initial configuration of rods\n#\n')
        seg_str = ""
        link_str = ""

        # Check for old data file to copy over
        if "data_file" in opts.params:
            data_path = Path(opts.params['data_file'])
            segments = read_sylinder_ascii_file(data_path)
            links = read_links_ascii_file(data_path)
            max_segment = max([int(seg.gid) for seg in segments])
            # Write segments to data file
            for seg in segments:
                seg_str += seg.to_string()
            for link in links:
                link_str += link.to_string()
            gen_gid = gen_id(int(max_segment) + 1)
        else:
            gen_gid = gen_id()

        start_pos = None
        for i, fil_params in enumerate(opts.params['filaments']):
            if (start_pos is not None and
                    not isinstance(fil_params.get('start_pos', None), list)):
                fil_params['start_pos'] = start_pos
            set_default_parameters(fil_params,
                                   opts.params['global_parameters'],
                                   check_for_none=True)

            print(f'### Making filament {i}')
            fil = FlexFilament(**fil_params)
            fil.make_flex_fil(gen_gid)
            seg_str += fil.get_fil_seg_str()
            link_str += fil.get_fil_link_str()

            if fil_params.get('join_next', False):
                start_pos = fil.get_next_fil_section_start(
                    fil_params['start_next_displace'])
            else:
                start_pos = None

        if fil_params.get('join_next', False):
            raise RuntimeError("Can't join next")

        f.write(seg_str)
        f.write(link_str)


##########################################
if __name__ == "__main__":
    opts = parse_args()
    gen_flexible_init(opts)
