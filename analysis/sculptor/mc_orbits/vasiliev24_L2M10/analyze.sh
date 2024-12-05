#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/vasiliev24_L2M10/
julia resample_lmc.jl $out_dir/trajlmc.txt
combine_outputs.py $out_dir/out
julia calc_peris.jl
ln -s $out_dir/initial.hdf5 .
