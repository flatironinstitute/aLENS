# _aLENS_ Utilities

This is a collection of scripts and modules to make creating, testing, and running _aLENS_ jobs easier.

## Initial data generator scripts for simulations

To access the full capabilities of _aLENS_, you will often need to create initial configurations for filaments, spheres, and crosslinking proteins. This can be done effectively by with scripts that algorithmically populate the space with proper number of filaments and proteins. Here are some example python with brief explanations of their purpose. See the doc strings in the files for more detailed use.

### Run generators

- `gen_flexible_init.py`: Generate flexibile filaments using a yaml file. (See flex_fil_param.yaml)
- `gen_long_filament_from_data`: Generate a long filament from the configurations of a short filament in time. Only works if the short filament was pinned at both ends so the end-to-end vector is constant in time.
- `gen_rigid_filament_protein_file.py`: Populate a filament with crosslinking proteins given a yaml file.

### Testing generators

- `gen_particle_cube.py`: Generates a cube made of 27 spheres. Used for unit testing.
- `scaling_sim_generator.py`: Create a series of identical simulations and sbatch scripts to test the timestep scaling as a function of nodes and cores. Simulations are based on the files in directory system is run in. (Currently only works for the Flatiron cluster architecture.)
