#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/vasiliev+21_lmc/
combine_outputs.py $out_dir/out
julia resample_lmc.jl $out_dir/trajlmc.txt
julia calc_peris.jl
ln -s $out_dir/initial.hdf5 .
