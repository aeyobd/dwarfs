#!/bin/bash

simdir="$DWARFS_SIMS/all_galaxies/M11"
rm simulation
rm combined.hdf5

ln -s $simdir simulation

combine_outputs.py simulation/out

julia ../extract_trajectory.jl .
