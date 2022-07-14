#!/usr/bin/env python

"""@package docstring
File: scaling_sim_generator.py
Author: Adam Lamson
Email: alamson@flatironinstitute.org
Description: Quick time scaling of very large simulations. Used to compare with
previous runs
"""

import os
import yaml
import shutil
import subprocess
from pathlib import Path

batch_script_str = """#!/bin/bash

###############################################################################
########################## SLURM Settings #####################################
#NOTE: if you request email, you will get one for each task!
# The stdout (-o) and stderr (-e) files will be written to
# THE SUBMISSION DIRECTORY with names that follow the given templates.
#SBATCH -o outrunLog.%j.%A.%a.%N.out
#SBATCH -e errrunLog.%j.%A.%a.%N.err

#SBATCH --job-name={0} 
#SBATCH --nodes={0}
#SBATCH --exclusive

# for one mpi rank per socket
#SBATCH --ntasks-per-socket=4
#SBATCH --cpus-per-task=16
#SBATCH --time 2:00:00

# Flatiron Institute specific
#SBATCH --constraint=rome
#SBATCH --partition=ccb
###############################################################################

source ./set_env.sh

export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "---------SLURM SETTINGS---------"
echo "Greetings from $SLURM_JOB_NUM_NODES!"
env | grep -i slurm
echo "------------------------------"

echo "---------OMP SETTINGS---------"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_DISPLAY_ENV=true
export OMP_NESTED=false
export OMP_MAX_ACTIVE_LEVELS=1
env | grep -i OMP_
echo "------------------------------"

echo "---------KMP SETTINGS---------"
env | grep -i KMP_
echo "------------------------------"

echo "---------MKL SETTINGS---------"
env | grep -i MKL_
echo "------------------------------"

echo "---------MPI SETTINGS---------"
which mpirun
which mpiexec
which srun
env | grep -i MPI_
echo "------------------------------"

# let mpirun auto detect mpi settings arranged by slurm
echo "running on hosts:"
echo $SLURM_JOB_NODELIST
mpirun --map-by ppr:$SLURM_NTASKS_PER_SOCKET:socket:PE=$SLURM_CPUS_PER_TASK \
       --bind-to core --display-map ./aLENS.X

date
"""


def main():
    # TODO Make modifications to the runconfig so that you can use it with any runconfig
    #      -Add arguments for number of time steps
    # with (Path.cwd() / 'RunConfig.yaml').open('r') as pf:
    #     run_params = yaml.safe_load(pf)
    nnodes = [2**i for i in range(5)]
    root_path = Path.cwd()
    file_list = [f for f in root_path.glob('*') if not f.is_dir()]
    for n in nnodes:
        # Create a directory that will submit a job using N nodes
        cur_dir = root_path / (root_path.name + f"_N{n}")
        cur_dir.mkdir(exist_ok=True)
        # Copy over all the files that are in the root directory
        for f in file_list:
            shutil.copy(f, str(cur_dir))
        job_script = cur_dir / 'submit_job.sh'
        # Create the
        with (job_script).open('w') as sf:
            sf.write(batch_script_str.format(n))

        # Submit jobs
        os.chdir(cur_dir)
        subprocess.call(['sbatch', str(job_script)])
        os.chdir(root_path)
        # TODO Generate python script to test time
        # TODO write quick analysis for time
        # TODO Collect data once runs finish and give a report


##########################################
if __name__ == "__main__":
    main()
