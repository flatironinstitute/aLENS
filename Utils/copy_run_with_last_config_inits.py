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

    # Create new init files from the last configuration
    for result_path in root_dir.glob('**/result(|.zip)'):
        sd_dir = result_path.parent
        print(sd_dir)
        args.root_dir = sd_dir
        make_last_step_new_init(args)

    try:
        # Copy over run folder  (skipping result.zip and result directories)
        subprocess.call(["rsync",
                         "-avz",
                         "--exclude=*.log*",
                         "--exclude=*result*",
                         "--exclude=*.out*",
                         "--exclude=*.h5*",
                         "--exclude=*TimeStepInfo.txt*",
                         "--exclude=*disbatch_files*",
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
