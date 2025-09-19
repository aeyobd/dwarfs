#!/bin/bash

source paths.sh

set -xe

ln -s $out_path simulation

echo out $out_path
ls $out_path
combine_outputs.py $out_path/out/
calc_centres.jl --r_max 1
combine_outputs.py $out_path/out -c centres.hdf5
mass_profiles.jl -k 1

# for painting stars
