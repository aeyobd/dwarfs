#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage $0 foldername"
    exit 1
fi

out_path=$DWARFS_ROOT/analysis/isolation/1e7_new/fiducial/

stellar_profiles_3d.jl $out_path $1/probabilities_stars.hdf5 -o $1/stellar_profiles_3d.hdf5 --scale ../halo-used.toml 

stellar_profiles.jl $out_path $1/probabilities_stars.hdf5 -o $1/stellar_profiles.hdf5 --scale ../halo-used.toml 
