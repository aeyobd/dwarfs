#!/bin/bash

out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/vasiliev24_L3M11/
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia resample_lmc.jl
julia calc_peris.jl
