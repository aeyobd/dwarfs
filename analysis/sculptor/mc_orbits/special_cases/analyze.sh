#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/special_cases/out
combine_outputs.py $out_dir
julia calc_peris.jl
cp $out_dir/../initial.hdf5 .