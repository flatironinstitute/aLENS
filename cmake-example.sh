#! /bin/bash

# remember to use Release and proper flags for production runs
# examples for gcc:
# for production run:
# -D CMAKE_BUILD_TYPE=Release \
# -D CMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" \
# for debug/development:
# -D CMAKE_BUILD_TYPE=Debug \
# -D CMAKE_CXX_FLAGS="-O0" \

# choose a proper install destination
# change SFTPATH to where all dependencies are installed

cmake \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_TEST=ON \
  -D CMAKE_CXX_FLAGS="-O3 -march=native -DNDEBUG" \
  -D CMAKE_INSTALL_PREFIX="/your/installation/location" \
  -D SFTPATH="/your/dependency/library/location" \
  ../
