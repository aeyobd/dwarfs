#!/bin/bash
out_dir=$DWARFS_SIMS/sculptor/mc_orbits/vasiliev24_L3M11_2x_special_cases/
combine_outputs.py $out_dir/out
julia resample_lmc.jl $out_dir/trajlmc.txt
julia calc_peris.jl
ln -s $out_dir/initial.hdf5 .
