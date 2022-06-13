![](Docs/source/images/aLENS_Logo_RGB.png)
# aLENS (a Living ENsemble Simulator)

The motivation, algorithm and examples are discussed in this paper:
[Towards the cellular-scale simulation of motor-driven cytoskeletal assemblies](https://elifesciences.org/articles/74160)

# Normal users (on a laptop or a single desktop)

Please use the precompiled docker image distributed through [DockerHub](https://hub.docker.com/r/wenyan4work/alens).
First install docker on you computer (windows users need to install WSL2 first).
Then:

```bash
docker pull wenyan4work/alens
```

This command pulls everything you need into your local computer.
Then follow the `README` file in the docker image to run `aLENS`.

This docker image also contains the full set of development environment (compiler and dependence libraries).
You can edit the code as you wish.

**Note 1**: this development toolchain in this docker image is based on [spack](https://github.com/spack/spack) virtual environment.

**Note 2**: Using docker images is convenient but has its limitations:

1. You can only use multi-thread parallelization on a single computer. MPI is not supported.
2. You should not generated data files within the docker image. Use the filesystem mapping feature of docker to write the generated files to the host filesystem outside the docker image.

**Note 3**: Theoretically it is possible to run the docker image using singularity and parallel it with mpi, controlled by slurm.
However it highly depends on your local toolchain and cluster set up.
Maybe it will work well for you but we are not able to provide general support for this use case.

# Executable input: `Config.yaml` and `Initial.dat`

The executable `aLENS.X` reads expects 4 input files:

- `RunConfig.yaml` specifies configuration for system and MTs.
- `ProteinConfig.yaml` specifies configuration and number for proteins.
- `TubuleInitial.dat` specifies initial configuration of MTs.
- `ProteinInitial.dat` specifies initial configuration of proteins.

You can go to the folder `Examples/MixMotorSliding` to see examples of these files.

The two `Config.yaml` files are necessary, but the two `Initial.dat` files are optional.  
There are three cases:

- Case 1. No `dat` file exists. In this case MTs and proteins will be generated according to the settings in `RunConfig.yaml` and `ProteinConfig.yaml`
- Case 2. `TubuleInitial.dat` file exists, but `ProteinInitial.dat` does not. 
  In this case MTs will be read from the `TubuleInitial.dat`, and the MT number & length settings in `RunConfig.yaml` will be ignored.
  Proteins will be generated according to the settings in `ProteinConfig.yaml`.
- Case 3. Both `TubuleInitial.dat` and `ProteinInitial.dat` files exits. 
  In this case MTs will be read from the `TubuleInitial.dat`, and the MT number & length settings in `RunConfig.yaml` will be ignored. 
  Proteins will be read from the file `ProteinInitial.dat`. `aLENS` will try to reconstruct the initial binding status according to `ProteinInitial.dat`. 
  If reconstruction fails for a certain protein, for example, if a protein is specified to bind some MT but the MT does not appear at the correct location, an error message will be printed out and this end (that an error appears) of this protein will be set to unbound and the program **continues**.

In general, Case 1 is good for initiating a simulation and Case 3 is good for continuing a simulation with saved data files. Case 2 is useful for some cases where the effect of protein on a given MT configuration is of interest.

# The installation folder structure

Once 'make install' finishes, you will get a folder structure like the following. 
Assume that your installation folder is located at ~/Run
```bash
~/Run/
├── aLENS.X              # the executable
├── gitversion.txt       # the git hashtag for the executable
├── result               # the folder where data is saved
│   ├── Clean.sh         # the script to remove all data
│   ├── PNG              
│   │   ├── MovieGen.sh  # the script to generate movie using png sequences
│   │   └── cleanpng.sh
│   ├── Result2PVD.py    # create meta-file for Paraview to load data
│   ├── compress.sh      # compress data file into 7z for archiving
│   └── uncompress.sh    # uncompress 7z files
└── scripts
    ├── DensityMT.py     # simple microtubule density calculator
    ├── jobsub_slurm.sh  # example of slurm job submission script
    └── run_ompi.sh      # example of mpi run on a single multi-core machine

3 directories, 11 files
```


# The minimal set of necessary files

In the minimal case, you need only three files and a folder to run the executable:
- one executable `aLENS.X`.
- two input configuration files `RunConfig.yaml` and `ProteinConfig.yaml`.
- one folder `result` for saved data files.

# Your first run

Use the provided configuration and initial to run your first simulation.

```bash
$ cp ./Examples/MixMotorSliding/* ~/Run/
$ cd ~/Run/
$ ./aLENS.X > ./log.txt
```

# Data organization

The program `aLENS.X` outputs to the folder `result`. 
`result` is at the same folder as `aLENS.X` itself.

It first writes a file `simBox.vtk`, which shows the simulation box as a simple rectangular box. For example:

```bash
$ cat ./result/simBox.vtk
# vtk DataFile Version 3.0
vtk file
ASCII
DATASET RECTILINEAR_GRID
DIMENSIONS 2 2 2
X_COORDINATES 2 float
0 20
Y_COORDINATES 2 float
0 1
Z_COORDINATES 2 float
0 1
CELL_DATA 1
POINT_DATA 8
```

Then the executable writes several different sequences of data files.

Two of them are human readable ascii files:

- `SylinderAscii_*.dat` are human readable data files of MTs (Sylinder = Spherocylinder). These files can be directly used as `TubultInitial.dat` by renaming.
- `ProteinAscii_*.dat` are human readable data files of proteins. These files can be directly used as `ProteinInitial.dat` by renaming.

Four of them are XML vtk format in base64 binary encoding. These are not human readable but can be conveniently loaded into `Paraview` for visualization or read by VTK (either python or cpp) for data processing.

- `Sylinder_*.pvtp` save data for MTs.
- `Protein_*.pvtp` save data for proteins.
- `ConBlock_*.pvtp` save data for collision and protein constraint blocks.

For explanation of these `pvtp` files, read the official guide of vtk file format: `https://kitware.github.io/vtk-examples/site/VTKFileFormats/#parallel-file-formats`.
In short, each `pvtp` file (parallel vtp) is a tiny index to a set of `vtp` files (serial vtp), which holds the actual data.
`aLENS.X` is written such that each MPI rank writes its own set of data to a unique `vtp` file.
Therefore the number of `vtp` files in each `pvtp` file index is equal to the number of MPI ranks.
The restriction is that the index `pvtp` file must appear in the same location as those `vtp` data files.

These sequeces of files are divided into different subfolders so each subfolder contains no more than a few thousand files.
This is due to the limitations of parallel file systems on some mpi clusters where saving a large number of files in a single directory destroys IO performance or even crashes the executable.

The data files are saved in different folders, but for postprocessing & visualization, in some cases there are some restrictions that all files of the same sequence must appear in the same folder otherwise the postprocessing or visualization program may fail to load the entire sequence.
The python script `Result2PVD.py` is used to handle this situation.
It creates `.pvd` files, which are indices to those `.pvtp` files and can be loaded by `Paraview`, so that all files belong to one sequence can be found in one place.
It is safe to run this script when `aLENS.X` is still running and writing data.

To run this script:

```bash
$ cd ./result/
$ python3 ./Result2PVD.py
$ ls ./*.pvd
./ConBlockpvtp.pvd  ./Proteinpvtp.pvd  ./Sylinderpvtp.pvd
```

# Further documentation of code internals

You need `doxygen` to generate html documents for the code internals.
Go to the root folder of `aLENS` and let `doxygen` to generate the documentation, according to the configuration file `Doxyfile`:

```bash
$ doxygen
```

Then you can open the file `doc/html/index.html` in your browser to read internal code document.
For example, if you have firefox installed:

```
$ firefox doc/html/index.html
```

# NORMAL USERS PLEASE STOP HERE

Do not read any further if you are a normal user.

# COMPILATION

For best performance you should utilize `MPI+OpenMP` parallelization in any cases where you need to run the job across multiple `numa` regions.
If you do not know what `numa` means, stop here and use the docker image instead.

Parallization with MPI requires compilation from source because every MPI library on a cluster must interact closely with the networking hardware and no generic precompiled binary executable file is able to achieve that.

You need a fully working MPI compiler, where the base compiler must fully support not only the `c++14` language features, but also the full `STL` defined in `c++14`.
Read this [webpage](https://en.cppreference.com/w/cpp/compiler_support#cpp14) about compiler information.
`gcc-5` could possibly work but we recommend `gcc-7` or higher.

You also need `cmake >= 3.10` to configure this project.
If you do not have it, download it from `https://cmake.org/download/`.

## Step 1, compile and install dependencies

You will have to install the following 4 packages:

1. `boost>=1.71`
2. `Eigen>=3.3`
3. `vtk>=9`
4. `Trilinos=12.18.1`

We provide reference cmake scripts under the folder `Dep` to compile those libraries in the correct mode.
For the performance of `Trilinos`, you should choose a proper linear algebra package **WITH MULTI-THREADING SUPPORT**, such as `OpenBLAS` or `MKL`.
Example cmake scripts have been provided for each case.
If you are using some other linear algebra libraries, modify the cmake script for `Trilinos` accordingly.

For convenience we provide two scripts under the folder `Dep` to automate this process.

```bash
$ cd Dep # go to Dep folder
$ python3 download_all.py # download source code of all dependence libraries
$ python3 compile_all.py # compile and install all dependence libraries
```

Before running `compile_all.py`, modify the `user switches` section of `compile_all.py` to set switches that match your toolchain.

## Step 2, compile and install `aLENS`

Once those dependence libraries have been installed, compile and install `aLENS` use the provided `cmake-example.sh` script.
Remember to set `SFTPATH` to the location where you installed those dependence libraries.

# `MPI+OpenMP` run

This is the actual running mode for a cluster.
However, due to the inconsistency of how `Intel mpi`, `openmpi`, and `mpich` handles multithread mapping, you may need different environment variable settings.
Here is a brief summary of things you may have to tune.

In the folder `Run/scripts` you can find `jobsub_slurm.sh` as an example of how to submit jobs to slurm with openmpi.

If you are not sure how you should setup things, consult your system administrator or play with [AMASK](https://github.com/TACC/amask) to see how different settings affect different threading mapping and binding modes.

## Environment Variable Settings

### OpenMP

- `OMP_NESTED=FALSE` **REQUIRED** and/or `OMP_MAX_ACTIVE_LEVELS=1` for new compilers.
- `OMP_NUM_THREADS=N` Change `N` to the number of cores per MPI rank. Hyperthreading may or may not be useful.
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

- `OPENBLAS_NUM_THREADS=N` Change `N` to thread-parallel tune performance. Read `OpenBLAS` document to get more information about this environment variable.
