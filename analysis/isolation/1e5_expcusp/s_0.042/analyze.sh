#!/bin/bash

source paths.sh
rm *.hdf5
rm *.csv

set -xe


ln -s $out_path/ simulation
combine_outputs.py simulation/out/
calc_centres.jl --r_max 10
combine_outputs.py simulation/out -c centres.hdf5
mass_profiles.jl -k 1
cp simulation/halo.toml .
