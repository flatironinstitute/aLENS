# Installing aLENS from scratch

<!-- Created By: Adam Lamson -->
<!-- Last Edited: August 16, 2022 -->

## General tasks to do beforehand if developing aLENS

- Get github account and collaborator status on aLENS and SimToolbox
- Install git and setup up token access in global git config
- (Optional) Get Intel license to install MKL on system

## Ubuntu (version 22)

**WARNING**: Make sure you have more than 2Gb of RAM and 30Gb of hard drive storage. There are work arounds for reduced memory and storage but these are not covered in this tutorial.

1. (Optional) If on virtual machine set up port forwarding to ssh into desktop
   [SSH into virtual machine in virtual box](https://bobcares.com/blog/virtualbox-ssh-nat/)
1. Upgrade to make sure you have latest `apt` version
   ```bash
   sudo apt update
   sudo apt upgrade
   ```
1. Install MKL on system (from this [blog](https://www.r-bloggers.com/2018/04/18-adding-intel-mkl-easily-via-a-simple-script/))

   ```bash
   wget -P /tmp [https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB](https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB)
   sudo apt-key add /tmp/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
   sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
   sudo apt update
   sudo apt install intel-mkl-64bit-2020.0-088
   ```

   `apt-key` command may throw a warning. This is fine, but for further reading, please see:
   [What commands (exactly) should replace the deprecated apt-key?](https://askubuntu.com/questions/1286545/what-commands-exactly-should-replace-the-deprecated-apt-key)

1. Tell computer where to find mkl by setting

   ```bash
   export MKLROOT=/opt/intel/mkl
   ```

   Linux usually installs the mkl package in `/opt` directory but check to make sure `intel/mkl` is there. Otherwise, go searching for it.
   Be sure to add `MKLROOT` variable to `.bashrc` file so that you don’t always have to type above command when you open a new terminal. Quick way to add this command to your path is with the command

   ```bash
   echo 'export MKLROOT=/opt/intel/mkl' >> ~/.bashrc
   ```

1. Install other necessary dependencies

   ```bash
   sudo apt install gcc g++ cmake openmpi-bin openmpi-doc libopenmpi-dev lbzip2
   ```

   **Make sure your `make` version <4.3 otherwise Trilinos12 [will not compile](https://github.com/UoB-HPC/BabelStream/issues/104).** This can be installed by first removing `make` and then downloading and installing the debian package

   ```bash
   sudo apt-get remove make
   wget http://ftp.de.debian.org/debian/pool/main/m/make-dfsg/make_4.2.1-1.2_amd64.deb
   sudo dpkg -i make_4.2.1-1.2_amd64.deb
   ```

1. Make directories to store program files
   ```bash
   mkdir -p ~/projects ~/local/aLENS
   ```
1. Create project structure and clone

   ```bash
   cd ~/projects
   git clone --recursive git@github.com:flatironinstitute/aLENS.git
   cd aLENS/Dep
   ```

1. Download depency softwares using the provided script

   ```bash
   python3 download_all.py
   ```

   This can take longer than you think. Just be patient.

1. Make necessary changes to `compile_all.py` file which includes changing installation directory to
   ```bash
   # your installation destination, use ABSOLUTE PATH
   ##From
   #os.environ["SFTPATH"] = os.environ['HOME'] + \
   #    '/envs/alens_env'
   ##To
   os.environ["SFTPATH"] = os.environ['HOME'] + \
       '/local/'
   ```
1. Compile all the dependencies
   ```bash
   python3 compile_all.py
   ```
   This will take a **_LONG_** while. Go take a break and come back.
1. (Now following the installation guide on [github landing page](https://github.com/flatironinstitute/aLENS/).) Go back to aLENS root directory and create a build directory

   ```bash
   cd ..
   mkdir build
   ```

1. Modify `cmake-example.sh` by replacing the variables `SFTPATH` and `CMAKE_INSTALL_PREFIX` definitions with
   ```bash
   -D CMAKE_INSTALL_PREFIX="${HOME}/local/aLENS" \
   -D SFTPATH="${HOME}/local" \
   ```
1. Go to build directory, run cmake script, and make executable

   ```bash
   cd build
   sh ../cmake-example.sh
   make -j <number_of_compiling_cores>
   ```

   **WARNING**: Leaving `<number_of_compiling_cores>` blank creates threads equal to the number of cores on your local machine. If you do this with a virtual machine, it will most likely crash during compilation.

   Compiling creates the executable file `aLENS.X` which you can run in any of the example simulations included in the `Examples` directory of aLENS.

## Flatiron cluster (2022-05-18)

This is much easier but only because all dependencies are already installed on our module system.

1. Make useful directories to store program files
   ```bash
   mkdir -p ~/projects ~/local/aLENS
   ```
1. Create project structure and clone
   ```bash
   cd ~/projects
   git clone --recursive git@github.com:flatironinstitute/aLENS.git
   cd aLENS
   ```
1. Load essential modules
   ```bash
   module purge
   module load modules gcc cmake gsl boost lib/fftw3 intel-mkl openmpi4 trilinos/12.18.1-mpi eigen vtk
   ```
1. Export environment variables
   ```bash
   export MKL_INTERFACE_LAYER=GNU,LP64
   export MKL_THREADING_LAYER=GNU
   export BOOST_ROOT=$BOOST_BASE
   export FFTWDIR=$FFTW3_BASE
   unset OMPI_CC
   unset OMPI_CXX
   export OMP_DISPLAY_ENV=true
   export OMP_MAX_ACTIVE_LEVELS=1
   ```
1. Modify `cmake-example.sh` by replacing SFTPATH and CMAKE_INSTALL_PREFIX definitions with
   ```bash
   -D CMAKE_INSTALL_PREFIX="${HOME}/local/aLENS" \
   -D SFTPATH="/cm/shared/sw/nix/store" \
   ```
1. Make and go to build directory, run cmake script, and make the executable
   ```bash
   mkdir build
   cd build
   sh ../cmake-example.sh
   make -j
   ```

This will produce an executable file `aLENS.X` that you can run in any of the simulation directory included in the `Examples` directory of aLENS.
