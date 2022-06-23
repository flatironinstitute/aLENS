# Installing aLENS from scratch

Created By: Adam Lamson
Last Edited: June 8, 2022 10:22 PM

## General tasks to do beforehand if developing

- Get github account and collaborator status on aLENS and SimToolbox
- Install git and setup up token access in global git config
- (Optional) Get Intel license to install MKL on system

## Ubuntu (version 22)

WARNING: Make sure you have more than 2Gb of RAM and 30Gb of hard drive storage. There are work arounds for reduced memory and storage but they are not covered in this worksheet currently.

- (Optional) If on virtual machine set up port forwarding to ssh into desktop
    
    [SSH into virtual machine in virtual box](https://www.notion.so/SSH-into-virtual-machine-in-virtual-box-a50f95a71020402ea9b05b65f0fdae45)
    
- Upgrade to make sure you have latest `apt` version
    
    ```bash
    sudo apt update
    sudo apt upgrade
    ```
    
- Install MKL on system (from this [blog](https://www.r-bloggers.com/2018/04/18-adding-intel-mkl-easily-via-a-simple-script/))
    
    ```bash
    wget -P /tmp [https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB](https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB)
    sudo apt-key add /tmp/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
    sudo apt update
    sudo apt install intel-mkl-64bit-2020.0-088
    ```
    
    `apt-key` command may throw warning. This is fine but for further reading
    
    [What commands (exactly) should replace the deprecated apt-key?](https://askubuntu.com/questions/1286545/what-commands-exactly-should-replace-the-deprecated-apt-key)
    
    Tell computer where to find mkl by setting
    
    ```bash
    export MKLROOT=/opt/intel/mkl
    ```
    
    Linux usually installs it in `/opt` directory but make sure `intel/mkl` is there otherwise go searching for it.
    
    Be sure to add `MKLROOT` variable to `.bashrc` file so that you don’t always have to type above command.  Quick way to add the path is
    
    ```bash
    echo 'export MKLROOT=/opt/intel/mkl' >> .bashrc
    ```
    
- Install other necessary dependencies
    
    ```bash
    sudo apt install gcc g++ cmake openmpi-bin openmpi-doc libopenmpi-dev lbzip2
    ```
    
    **Make sure your `make` version <4.3 otherwise Trilinos12 [will not compile](https://github.com/UoB-HPC/BabelStream/issues/104).** This can be installed by first removing `make` and then downloading and installing the debian package 
    
    ```bash
    sudo apt-get remove make
    wget http://ftp.de.debian.org/debian/pool/main/m/make-dfsg/make_4.2.1-1.2_amd64.deb
    sudo dpkg -i make_4.2.1-1.2_amd64.deb
    ```
    
- Make useful directories to store program files
    
    ```bash
    mkdir -p ~/projects ~/local/aLENS
    ```
    
- Create project structure and clone
    
    ```bash
    cd ~/projects
    git clone --recursive git@github.com:flatironinstitute/aLENS.git
    cd aLENS/Dep
    ```
    
- Download files for dependencies to install
    
    ```bash
    python3 download_all.py
    ```
    
    This can take longer than you think. Just be patient.
    
- Make necessary changes to `compile_all.py` file which includes changing installation directory to
    
    ```bash
    # your installation destination, use ABSOLUTE PATH
    ##From
    #os.environ["SFTPATH"] = os.environ['HOME'] + \
    #    '/envs/alens_env'
    ##To
    os.environ["SFTPATH"] = os.environ['HOME'] + \
        '/local/'
    ```
    
- Compile all the dependencies
    
    ```bash
    python3 compile_all.py
    ```
    
    This will take a ***LONG*** while. Go take a break and come back.
    
- Follow remaining installation guide on github landing page. First go back to aLENS directory and create a build directory
    
    ```bash
    cd ..
    mkdir build
    ```
    
- Modify `[cmake-example.sh](http://cmake-example.sh)` by replacing SFTPATH and CMAKE_INSTALL_PREFIX definitions with
    
    ```bash
    -D CMAKE_INSTALL_PREFIX="${HOME}/local/aLENS" \
    -D SFTPATH="${HOME}/local" \
    ```
    
- Go to build directory, run cmake script, and make executable
    
    ```bash
    cd build
    sh ../cmake-example.sh
    make -j <number_of_compiling_cores>
    ```
    

Note: leaving `<number_of_compiling_cores>` blank creates threads equal to the number of cores on your local machine. If you do this with a virtual machine, it will most likely crash during compilation. 

Compiling creates the executable file `aLENS.X` which you can run in any of the example simulations included in the `Examples` directory of aLENS.

## Mac installation guide

TODO… This is a pain

## Flatiron cluster (2022-05-18)

This is much easier but only because of all dependencies are installed on our module system

- Make useful directories to store program files
    
    ```bash
    mkdir -p ~/projects ~/local/aLENS
    ```
    
- Create project structure and clone
    
    ```bash
    cd ~/projects
    git clone --recursive git@github.com:flatironinstitute/aLENS.git
    cd aLENS
    ```
    
- Load essential modules
    
    ```bash
    module purge
    module load modules gcc cmake gsl boost lib/fftw3 intel-mkl openmpi4 trilinos/12.18.1-mpi eigen vtk1
    ```
    
- Export environment variables
    
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
    
- Modify `[cmake-example.sh](http://cmake-example.sh)` by replacing SFTPATH and CMAKE_INSTALL_PREFIX definitions with
    
    ```bash
    -D CMAKE_INSTALL_PREFIX="${HOME}/local/aLENS" \
    -D SFTPATH="/cm/shared/sw/nix/store" \
    ```
    
- Make and go to build directory, run cmake script, and make the executable
    
    ```bash
    mkdir build
    cd build
    sh ../cmake-example.sh
    make -j
    ```

This will produce an executable file `aLENS.X` which you can run in any of the example simulations included in the `Examples` directory of aLENS.
