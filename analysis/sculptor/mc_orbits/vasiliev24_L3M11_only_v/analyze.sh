#!/bin/bash

out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/vasiliev24_L3M11_only_v/
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia resample_lmc.jl
julia calc_peris.jl
