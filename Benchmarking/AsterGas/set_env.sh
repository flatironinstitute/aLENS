#!/usr/bin/env bash

################################################################################
#Name:
#Description:
#Input:
#Output:
################################################################################

echo $1

module purge

# unset redundant envs
#unset $(compgen -v | grep "MKL_$")
#unset $(compgen -v | grep "OMP_$")
unset MKLROOT
unset MKL_INTERFACE_LAYER
unset MKL_THREADING_LAYER
unset MKL_INCLUDE_DIRS
unset MKL_LIB_DIRS
unset FFLAGS
unset CFLAGS
unset CXXFLAGS

# toolchain
module load modules
module load slurm
module load cmake
module load gcc

echo "loading modules."
module load intel-mkl
module load openmpi4

# library
module load gsl
module load boost
module load eigen
module load vtk
module load lib/fftw3
module load trilinos/mpi-12.18.1

# Reset necessary environment variables
export MKL_INTERFACE_LAYER=GNU,LP64
export MKL_THREADING_LAYER=GNU
export BOOST_ROOT=$BOOST_BASE
export FFTWDIR=$FFTW3_BASE
unset OMPI_CC
unset OMPI_CXX
export OMP_DISPLAY_ENV=true
export OMP_MAX_ACTIVE_LEVELS=1

# set anaconda environment if you wish
#conda activate alens
