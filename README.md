# aLENS (active Living ENsemble Simulator)
The motivation, algorithm and examples are discussed in this paper: 
[aLENS: towards the cellular-scale simulation of motor-driven cytoskeletal assemblies](https://arxiv.org/abs/2109.08206)

# Introduction

This file will guide you through the preparation, compilation, running, and postprocessing of using aLENS. 

To use this software, the minimum requirement is a mpi c++ compiler fully supporting `c++14` and `openmp-4.0`, for example, `mpicxx` that calls `gcc>=7`.

# Clone the repo
First, clone this repo:
```bash
git clone https://github.com/wenyan4work/aLENS.git
```

This repo relies on two submodules. You need to initialize them after cloning:
```bash
cd aLENS
git submodule init
git submodule update
```

# Install dependencies
## Step 1. dependencies of SimToolbox
An easier way is to use the automated compilation and installation scripts hosted at `https://github.com/wenyan4work/Environment`.

## Step 2. dependencies of KMC 
KMC relies on the Gauss-Kronrod quadrature integrator in Boost.math. Make sure you have `boost>=1.70` installed.
Only the boost header is necessary. 


# Compiling
Use 'cmake' to compile things. 
If `make` finishes successfully, now you have a usable executable called `aLENS.X`. 


# Executable input: `Config.yaml` and `Initial.dat`
The executable `aLENS.X` reads expects 4 input files:
- `RunConfig.yaml` specifies configuration for system and MTs.
- `ProteinConfig.yaml` specifies configuration and number for proteins.
- `TubuleInitial.dat` specifies initial configuration of MTs.
- `ProteinInitial.dat` specifies initial configuration of proteins.

You can go to the folder `InitialConfigExample/SimplePair` to see examples of these files.

The two `Config.yaml` files are necessary, but the two `Initial.dat` files are optional.  
There are three cases:
- Case 1. No `dat` file exists. In this case MTs and proteins will be generated according to the settings in `RunConfig.yaml` and `ProteinConfig.yaml`
- Case 2. `TubuleInitial.dat` file exists, but `ProteinInitial.dat` does not. In this case MTs will be read from the `TubuleInitial.dat`, and the MT number & length settings in `RunConfig.yaml` will be ignored. Proteins will be generated according to the settings in `ProteinConfig.yaml`.
- Case 3. Both `TubuleInitial.dat` and `ProteinInitial.dat` files exits. In this case MTs will be read from the `TubuleInitial.dat`, and the MT number & length settings in `RunConfig.yaml` will be ignored. Proteins will be read from the file `ProteinInitial.dat`. `aLENS` will try to reconstruct the initial binding status according to `ProteinInitial.dat`. If reconstruction fails for a certain protein, for example, if a protein is specified to bind some MT but the MT does not appear at the correct location, an error message will be printed out and this end (that an error appears) of this protein will be set to unbound and the program continues.

In general, Case 1 is good for initiating a simulation and Case 3 is good for continuing a simulation with saved data files. Case 2 is useful for some cases where the effect of protein on a given MT configuration is of interest.

# The minimum set of necessary files
In the minimu case, you need only three files and a folder to run the executable:
- one executable `aLENS.X`.
- two input configuration files `RunConfig.yaml` and `ProteinConfig.yaml`.
- one folder `result` for saved data files.

# Your first run
Use the provided example to run your first simulation:
```bash
wyan$ cp ./InitialConfigExample/PairBinding/* ./
wyan$ ./aLENS.X > ./outrun.log
```

# Data organization
The program `aLENS.X` outputs to the folder `result`. `result` is at the same folder as `aLENS.X` itself.

It first writes a file `simBox.vtk`, which shows the simulation box as a simple rectangular box. For example:
```bash 
wyan$ cat ./result/simBox.vtk 
# vtk DataFile Version 3.0
vtk file
ASCII
DATASET RECTILINEAR_GRID
DIMENSIONS 2 2 2
X_COORDINATES 2 float
0 10
Y_COORDINATES 2 float
0 10
Z_COORDINATES 2 float
0 10
CELL_DATA 1
POINT_DATA 8
```
Then the executable writes 6 different sequences of data files. 

Two of them are human readable ascii files:
- `SylinderAscii_*.dat` are human readable data files of MTs (Sylinder = Spherocylinder). These files can be directly used as `TubultInitial.dat` by renaming.
- `ProteinAscii_*.dat` are human readable data files of proteins. These files can be directly used as `ProteinInitial.dat` by renaming.
  
Four of them are XML vtk format in base64 binary encoding. These are not human readable but can be conveniently loaded into `Paraview` for visualization or read by VTK (either python or cpp) for data processing.
- `Sylinder_*.pvtp` save data for MTs.
- `Protein_*.pvtp` save data for proteins.
- `ConBlock_*.pvtp` save data for collision and protein constraint blocks.

For explanation of these `pvtp` files, read the official guide of vtk file format: `https://lorensen.github.io/VTKExamples/site/VTKFileFormats/#parallel-file-formats`. 
In short, each `pvtp` file (parallel vtp) is a tiny index to a set of `vtp` files (serial vtp), which holds the actual data.
`aLENS.X` is written such that each MPI rank writes its own set of data to a unique `vtp` file. 
Therefore the number of `vtp` files in each `pvtp` file index is equal to the number of MPI ranks. 
The restriction is that the index `pvtp` file must appear in the same location as those `vtp` data files.

These sequeces of files are divided into different subfolders so each subfolder contains no more than roughly 3000 files. 
This is due to the limitations of parallel file systems on some mpi clusters where saving a large number of files in a single directory destroys IO performance or even crashes the executable.

The data files are saved in different folders, but for postprocessing & visualization, in some cases there are some restrictions that all files of the same sequence must appear in the same folder otherwise the postprocessing or visualization program may fail to load the entire sequence. 
The python script `Result2PVD.py` is used to handle this situation. 
It creates `.pvd` files, which are indices to those `.pvtp` files and can be loaded by `Paraview`, so that all files belong to one sequence can be found in one place. 
It is safe to run this script when `aLENS.X` is still running and writing data.

To run this script:
```bash
wyan$ cd ./result/
wyan$ python3 ./Result2PVD.py 
wyan$ ls ./*.pvd
./ConBlockpvtp.pvd  ./Proteinpvtp.pvd  ./Sylinderpvtp.pvd
```

# `MPI+OpenMP` run

This is the actual running mode for a cluster.
However, due to the inconsistency of how `Intel mpi`, `openmpi`, and `mpich` handles multithread mapping, you may need different environment variable settings.
Here is a brief summary of things you may have to tune.

In the folder `scripts` you can find `jobsub.slurm` as an example of how to submit jobs to slurm with openmpi.

If you are not sure how you should setup things, consult your system administrator or play with [AMASK](https://github.com/TACC/amask) to see how different settings affect different threading mapping and binding modes. 

## Environment Variable Settings
### OpenMP
- `OMP_NESTED=FALSE` **REQUIRED** and/or `OMP_MAX_ACTIVE_LEVELS=1` for new compilers.
- `OMP_NUM_THREADS=N` Change `N` to the number of cores. Hyperthreading may or may not be useful.
- `OMP_DISPLAY_ENV=VERBOSE` **Recommended**, helpful for debugging environment variables
- `OMP_PROC_BIND=spread` **Recommended by Kokkos**, you should check if it works for you
- `OMP_PLACES=threads` **Recommended by Kokkos**, you should check if it works for you


### If using Intel MKL:
- `MKL_THREADING_LAYER=INTEL` or `MKL_THREADING_LAYER=GNU` depending on if the compiler is `icpc` or `g++`
- `MKL_INTERFACE_LAYER=LP64` **Never** use `ILP64`

### If using Intel MPI (optional):
- `I_MPI_ADJUST_ALLTOALLV=1`
- `I_MPI_ADJUST_ALLTOALL=1`

### If using OpenBLAS (optional):
- `OPENBLAS_NUM_THREADS=N` Change `N` to tune multithreading. .

The priorities are `OPENBLAS_NUM_THREADS` > `GOTO_NUM_THREADS` > `OMP_NUM_THREADS`.
If you compile OpenBLAS with `USE_OPENMP=1`, you should set the `OMP_NUM_THREADS` environment variable; OpenBLAS ignores `OPENBLAS_NUM_THREADS` and `GOTO_NUM_THREADS` when compiled with `USE_OPENMP=1`.