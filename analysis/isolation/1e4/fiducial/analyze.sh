#!/bin/bash

source paths.sh

rm combined.hdf5 simulation centres.hdf5 profiles.hdf5 
set -xe

ln -s $out_path simulation
combine_outputs.py $out_path/out/
calc_centres.jl --r_max 10
combine_outputs.py $out_path/out -c centres.hdf5
mass_profiles.jl -k 1
