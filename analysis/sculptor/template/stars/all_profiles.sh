#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage $0 foldername"
    exit 1
fi
source ../paths.sh


stellar_profiles_3d.jl $iso_path $1/probabilities_stars.hdf5 -o $1/stellar_profiles_3d.hdf5 --scale ../halo-used.toml 

stellar_profiles.jl $iso_path $1/probabilities_stars.hdf5 -o $1/stellar_profiles.hdf5 --scale ../halo-used.toml 
