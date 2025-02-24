# Simulations


This directory contains the minimum amount of code to run each simulation and is where the raw simulation outputs should be stored.

- `isolation`. Because gravity is scale free, we only need a single suite of isolation simulations
- `sculptor` Gadget simulations of sculptor galaxies in a milky way potential.
- `ursa_minor` Gadget simulations of Ursa Minor galaxies in a milky way potential.


Each of these directories contains sub-directories for each halo (specified by resolution and properties). These directories then each contain subdirectories for different simulation parameters (e.g. orbits, potentials, gadget force accuracy, etc.). 

Inside each simulation directory, there should be a minimum of the following
- `run.sh` a script that runs the simulation
- `param.txt` the Gadget4 parameter file. These should all be a minimum change to template_param.txt in this folder
- `make_init.sh` creates the initial conditions. Note that these initial conditions may depend on another simulation
- `agama.ini` the agama potential file for the simulation

For a completed run, there will additionally be the following in a directory
- `out` contains the Gadget-4 simulation outputs, including status files on memory, cpu, etc., and every output snapshot
- `initial.hdf5` the initial conditions used in the simulation.
- `log.out` is the bash output of Gadget.
- `<#>.log` are the slurm bash logs for each run.


## Running simulations


