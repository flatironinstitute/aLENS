#! /bin/bash

# remember to use Release and proper flags for production runs
# examples for gcc:
# Production
cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" \
  -D CMAKE_INSTALL_PREFIX="/your/installation/location" \
  -D SFTPATH="/your/dependency/library/location" \
../

## Debug
# cmake \
#  -D CMAKE_CXX_COMPILER=mpicxx \
#  -D CMAKE_C_COMPILER=mpicc \
#  -D CMAKE_BUILD_TYPE=Debug \
#  -D ENABLE_TEST=ON \
#  -D CMAKE_CXX_FLAGS="-O0 -march=native -DDEBUG" \
#  -D CMAKE_INSTALL_PREFIX="/your/installation/location" \
#  -D SFTPATH="/your/dependency/library/location" \
# ../
