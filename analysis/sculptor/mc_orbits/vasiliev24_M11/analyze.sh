#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/vasiliev24_M11

bash clean.sh
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia ../calc_peris.jl
ln -s simulation/initial.hdf5 .
