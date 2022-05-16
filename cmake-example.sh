#! /bin/bash

# remember to use Release and proper flags for production runs
# examples for gcc:
# Production
cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS="-O3 -march=broadwell -DNDEBUG" \
  -D CMAKE_INSTALL_PREFIX="${HOME}/Run/AMSOS" \
  -D SFTPATH="${HOME}/local/" \
  -D Eigen3_DIR="${SFTPATH}/share/eigen3/cmake" \
../

# cmake \
#  -D CMAKE_CXX_COMPILER=mpicxx \
#  -D CMAKE_C_COMPILER=mpicc \
#  -D CMAKE_BUILD_TYPE=Debug \
#  -D ENABLE_TEST=ON \
#  -D CMAKE_CXX_FLAGS="-O0 -march=broadwell -DDEBUG" \
#  -D CMAKE_INSTALL_PREFIX="${HOME}/Run/AMSOS" \
#  -D SFTPATH="${HOME}/local/" \
#  -D Eigen3_DIR="${SFTPATH}/share/eigen3/cmake" \
# ../
