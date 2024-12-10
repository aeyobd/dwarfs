#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/EP20/out/
combine_outputs.py $out_dir
julia calc_peris.jl
rm initial.hdf5
ln -s $out_dir/../initial.hdf5 .
