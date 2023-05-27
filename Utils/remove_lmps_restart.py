#!/usr/bin/env python

"""@package docstring
File: copy_run_with_last_config_inits.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

from pathlib import Path
import re
# from make_last_step_new_init import make_last_step_new_init, get_args
import subprocess
from shutil import copy

# Find all directories with result.zip or result directory


def restart_number(f):
    """TODO: Docstring for restart_number.

    @param arg1 TODO
    @return: TODO

    """
    pattern = re.compile(r'.*restart\.([0-9]*)')
    if isinstance(f, str):
        return int(pattern.findall(f)[-1])

    return int(pattern.findall(f.name)[-1])


def main():
    root_dir = Path.cwd()

    # Try to find the correct result paths
    result_paths = list(root_dir.glob(r'**/result.zip'))
    if len(result_paths) == 0:
        result_paths = list(root_dir.glob(r'**/result/'))
    if len(result_paths) == 0:
        raise FileNotFoundError(
            'No result files/folders were found. May need to switch from "result.zip" to "result" to find the correct files.')
        return

    # Create new init files from the last configuration
    for result_path in result_paths:
        sd_dir = result_path.parent
        print(result_path)
        restart_paths = sorted(
            list(
                sd_dir.glob(r'*restart.[0-9]*')),
            key=restart_number)
        copy(restart_paths[-1], sd_dir / 'filament.restart.A')
        copy(restart_paths[-2], sd_dir / 'filament.restart.B')
        for restart_p in restart_paths:
            restart_p.unlink()


##########################################
if __name__ == "__main__":
    main()
