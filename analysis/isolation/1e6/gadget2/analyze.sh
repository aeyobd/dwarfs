#!/bin/bash

source paths.sh

set -xe

combine_outputs.py $out_path/out/
calc_centres.jl --r_max 1
combine_outputs.py $out_path/out -c centres.hdf5
mass_profiles.jl -k 1

cp $out_path/halo.toml .

get_energies.jl combined.hdf5/$idx_stars energies.hdf5
calc_dist_ana.jl energies.hdf5 halo.toml distribution_function.hdf5

