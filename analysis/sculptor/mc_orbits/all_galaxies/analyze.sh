#!/bin/bash

simdir="$DWARFS_ROOT/simulations/sculptor/mc_orbits/all_galaxies"
combine_outputs.py $simdir/out

resample_orbit.jl $simdir/trajlmc.txt -t combined.hdf5 -o orbit_lmc.csv --time-scale -1
