Post-processing and analysis
============================

.. _analysis:

Installation of aLENS_analysis
------------------------------

#. Clone alens_analysis on local machine::

    git clone --recursive https://github.com/flatironinstitute/aLENS_analysis.git
    cd aLENS_analysis
#. Create a virtual environment::

    conda create env -f environment.yml
    #or
    conda create env -n alens 
    conda install numpy h5py scipy matplotlib vtk pyyaml numba
    #or
    python3 -m venv /path/to/new/virtual/environment
    source my_venv/bin/activate
    pip install numpy h5py scipy matplotlib vtk pyyaml numba
#. Then install the alens_analysis::

    pip install -e .

Creating HDF5 data file from raw data
-------------------------------------

.. note::
    Need to go into this in more detail





