#!/bin/bash

###############################################################################
########################## SLURM Settings #####################################
#NOTE: if you request email, you will get one for each task!
# The stdout (-o) and stderr (-e) files will be written to
# THE SUBMISSION DIRECTORY with names that follow the given templates.
#SBATCH -o outrunLog.%j.%A.%a.%N.out
#SBATCH -e errrunLog.%j.%A.%a.%N.err

#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --export=ALL

# type of node
#SBATCH --constraint=skylake
# configuration per socket
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=20
###############################################################################

# remove old err and out files
find ./ -type f -name 'outrunLog*' -a ! -name "outrunLog.$SLURM_JOBID.*" -delete
find ./ -type f -name 'errrunLog*' -a ! -name "errrunLog.$SLURM_JOBID.*" -delete

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
echo "running on hosts:"
echo $SLURM_JOB_NODELIST
echo "------------------------------"

# let mpirun auto detect mpi settings arranged by slurm
# works for openmpi
mpirun --map-by ppr:$SLURM_NTASKS_PER_SOCKET:socket:PE=$SLURM_CPUS_PER_TASK \
       --bind-to core --display-map aLENS.X

date
