#!/bin/bash

source paths.sh
bash clean.sh
rm simulation

set -xe

ln -s $out_path simulation

combine_outputs.py simulation/out/
calc_centres.jl --r_max 1.0 > centres.log 2>&1
combine_outputs.py simulation/out -c centres.hdf5
mass_profiles.jl -k 1 > profiles.log 2>&1

cp simulation/orbit.csv .
