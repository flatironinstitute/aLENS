#!/usr/bin/env python

"""@package docstring
File: copy_run_with_last_config_inits.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description:
"""

from pathlib import Path
from make_last_step_new_init import make_last_step_new_init, get_args
import subprocess

# Find all directories with result.zip or result directory


def main():
    root_dir = Path.cwd()
    targ_dir = root_dir.parent / 'test_dir'

    args = get_args()

    args.step_info = True
    args.string_append = 'New'

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
        args.root_dir = sd_dir
        make_last_step_new_init(args)
        found_result_files = True

    try:
        # Copy over run folder  (skipping result.zip and result directories)
        subprocess.call(["rsync",
                         "-avz",
                         "--exclude=*.log*",
                         "--exclude=*result*",
                         "--exclude=*analysis*",
                         "--exclude=*.out*",
                         "--exclude=*.err*",
                         "--exclude=*.h5*",
                         "--exclude=*TimeStepInfo.txt*",
                         "--exclude=*disbatch_files*",
                         "--exclude=*TubuleInitial.dat*",
                         str(root_dir) + '/',
                         str(targ_dir)])

        # Move new initial files in target run directory to correct names
        for new_init_path in targ_dir.glob('**/TubuleInitialNew.dat'):
            sd_dir = new_init_path.parent
            new_init_path.rename(sd_dir / 'TubuleInitial.dat')
        for new_init_path in targ_dir.glob('**/ProteinInitialNew.dat'):
            sd_dir = new_init_path.parent
            new_init_path.rename(sd_dir / 'ProteinInitial.dat')

    except BaseException:
        raise

    finally:
        # Cleanup init paths left in root directories
        for new_init_path in root_dir.glob('**/*New.dat'):
            new_init_path.unlink()


##########################################
if __name__ == "__main__":
    main()
