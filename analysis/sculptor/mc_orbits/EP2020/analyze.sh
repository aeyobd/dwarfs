#!/bin/bash
bash clean.sh

out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/EP2020/
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia ../calc_peris.jl
