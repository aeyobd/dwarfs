#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/EP20_special_cases/
rm simulation
ln -s $out_dir simulation

combine_outputs.py simulation/out
julia calc_peris.jl
