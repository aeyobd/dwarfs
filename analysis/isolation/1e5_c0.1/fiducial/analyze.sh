#!/bin/bash

source paths.sh


bash clean

set -xe
ln -s $out_path simulation

combine_outputs.py simulation/out/
calc_centres.jl --r_max 1
combine_outputs.py simulation/out -c centres.hdf5
mass_profiles.jl -k $profile_skip

if [ "$copy_halo" = true ] ; then
    cp $out_path/halo.toml .
fi
