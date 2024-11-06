#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/all_studies/out
# combine_outputs.py $out_dir
echo calculating pericentres
julia calc_peris.jl

echo copying initial
cp $out_dir/../initial.hdf5 .
