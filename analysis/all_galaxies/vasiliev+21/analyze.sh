#!/bin/bash

simdir="$DWARFS_ROOT/simulations/all_galaxies/vasiliev+21"
ln -s $simdir simulation

combine_outputs.py simulation/out

resample_orbit.jl simulation/trajlmc.txt -t combined.hdf5 -o orbit_lmc.csv --time-scale -1
