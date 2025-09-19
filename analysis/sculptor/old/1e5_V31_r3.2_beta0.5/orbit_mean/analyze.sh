#!/bin/bash

source paths.sh
ln -s $out_path simulation

set -xe


combine_outputs.py $out_path/out/
calc_centres.jl --r_max 1.0
combine_outputs.py $out_path/out -c centres.hdf5
mass_profiles.jl -k 1
