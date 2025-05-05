#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/vasiliev24_L3M11_2x/

bash clean.sh
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia ../resample_lmc.jl simulation/trajlmc.txt
julia calc_peris_lmc.jl
ln -s simulation/initial.hdf5 .
