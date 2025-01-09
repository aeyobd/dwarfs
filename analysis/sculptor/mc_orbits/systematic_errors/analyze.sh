#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/vasiliev24_L2M11_loose/
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia calc_peris.jl
