#!/bin/bash

simdir="$DWARFS_SIMS/all_galaxies/L3M11_hhr"
rm simulation
rm combined.hdf5

ln -s $simdir simulation

combine_outputs.py simulation/out

julia ../extract_trajectory.jl .
