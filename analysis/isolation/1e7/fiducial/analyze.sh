#!/bin/bash

source paths.sh

rm centres.hdf5 combined.hdf5 profiles.hdf5 simulation halo.toml
ln -s $out_path simulation

set -xe

combine_outputs.py $out_path/out/
calc_centres.jl --r_max 1 2>&1 | tee centres.log
combine_outputs.py $out_path/out -c centres.hdf5
mass_profiles.jl -k 1 2>&1 | tee profiles.log

# cp $out_path/halo.toml .
