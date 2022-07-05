# First run of aLENS on a cluster (for dummies)

Created date: June 23, 2022 2:48 PM

This document is based on Adam Lamson’s notes and instructions. 

**Goal:** run one of the examples of aLENS on Rusty (FI cluster). 

# Connecting to a cluster

You need an account in FIDO. Once that is set, connect to any wifi (no need for FI wifi) and type into the terminal

```bash
ssh -p 61022 uursic@gateway.flatironinstitute.org
```

- It asks for a Google authentication code. That can be set up with a Google authenticator on your mobile device.
- then enter your FIDO password

This should connect you to the gateway. To connect to Rusty, type in  

```bash
ssh rusty
```

This brings you to your home directory. In my case, that is

```bash
/mnt/home/uursic
```

**NOTE:** don’t save (large) data into your home directory. The directory for saving data is

```bash
/mnt/ceph/users/uursic
```

## Setting up git and GitHub

**NOTE:** The default git version is 1.8.3.1, but there is a newer version available, just type in

```bash
module load git
```

You need to have a GitHub account ready. 

Setup your user name and user email address (use the email address you use for GItHub):

 

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@what.ever"
```

To connect your git on the cluster with your GitHub account, you need to get an SSH key. 

The next command generates a public and a private SSH key:

```bash
ssh-keygen -t ed25519
```

Press enter three times in order to avoid using a passphrase every time you want to log in. 

It tells you that your public and private keys have been generated and where you can access them. Type 

```bash
cat ~/.ssh/id_ed25519.pub
```

to access the public key. 

**NOTE:** never ever give your private key to anyone.

Go to your GitHub account online. Go to **settings > SSH and GPG keys** and click on the New SSH key. Copy the public key from terminal to your online GitHub account and give it a meaningful name. 

## Getting the aLENS repo from GitHub

In order to get everything into your local repository (on the cluster), make a new directory for aLENS and clone aLENS repository from GitHub:

```bash
mkdir -p ~/projects ~/local/aLENS
cd ~/projects
git clone --recursive git@github.com:flatironinstitute/aLENS.git
cd aLENS
```

comment: If `--recursive` is specified, this command will recurse into the registered submodules, and update any nested submodules within. 

Now is the time to load essential modules: 

```bash
module purge
module load modules gcc cmake gsl boost lib/fftw3 intel-mkl openmpi4 trilinos/12.18.1-mpi eigen vtk
```

comment: module purge unloads all the modules, except the sticky ones. If you call the command `module list`, you can get a list of all the loaded modules. The ones with (S) are sticky.

**NOTE:** loading vtk1 didn’t work for me, but loading vtk did (a difference between this and Adam’s notes)

The next step is to export environment variables

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

A way to do this is to make an .sh file, 

```bash
nano set_env.sh
```

paste in the text above and run it with an sh command:

```bash
sh set_env.sh
```

Modify `[cmake-example.sh](http://cmake-example.sh)` by replacing SFTPATH and CMAKE_INSTALL_PREFIX definitions with

```bash
-D CMAKE_INSTALL_PREFIX="${HOME}/local/aLENS" \
-D SFTPATH="/cm/shared/sw/nix/store" \
```

You can do this by opening the file with nano (or any other text editor):

```bash
nano cmake-example.sh
```

Make and go to build directory, run cmake script, and make the executable

```bash
mkdir build
cd build
sh ../cmake-example.sh
make -j
```

comment: -j specifies the number of jobs (commands) to run simultaneously.  If there is more than one -j option, the last one is effective.  If the -j option is given without an argument, make will not limit the number of jobs that can run simultaneously.

This will produce an executable file `aLENS.X` which you can run in any of the example simulations included in the `Examples` directory of aLENS.

## Running one of the examples

this chapter needs to be reviewed

Choose a directory where you want to run the code. In my case:

```bash
cd ~/local/aLENS
```

Now we need to copy a few things into this directory:

```bash
cp -r ~/projects/aLENS/Examples/MixMotorSliding .
cd MixMotorSliding
```

comment: The -r option just means that source directories will be copied as well as normal files.

and

```bash
cp -r ~/projects/aLENS/Run/* .
```

Earlier, we built the executable **aLENS.X** file. We will need it here as well. Let’s copy it into the directory

```bash
cp -r ~/projects/aLENS/build/aLENS.X ./aLENS.X
```

We are ready to run the executable now:

```bash
./aLENS.X
# or to control the number of cores used
OMP_NUM_THREADS=<number_of_cores> ./aLENS.
```

comment: I didn’t specify the number of cores for my first run and it took a long time to do anything. —> I need to figure this out. 

---

**General comment:** The folder structure here could have been much better. I will think of better ways to organize the files.