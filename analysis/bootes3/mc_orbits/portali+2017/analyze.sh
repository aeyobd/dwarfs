#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/bootes3/mc_orbits/portali+2017/
rm simulation
ln -s $out_dir simulation

combine_outputs.py $out_dir/out
julia calc_peris.jl
