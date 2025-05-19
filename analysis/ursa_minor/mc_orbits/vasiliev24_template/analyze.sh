#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/vasiliev24_FILL_IN

bash clean.sh
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia ../resample_lmc.jl simulation/trajlmc.txt
julia ../calc_peris_lmc.jl
ln -s simulation/initial.hdf5 .
