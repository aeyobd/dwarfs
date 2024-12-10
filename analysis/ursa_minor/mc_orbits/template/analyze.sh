#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/FILL_IN
combine_outputs.py $out_dir
julia calc_peris.jl
ln -s $out_dir/../initial.hdf5 .
