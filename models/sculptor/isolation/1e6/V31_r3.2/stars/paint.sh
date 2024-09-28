#!/bin/bash


if [ $# -ne 1 ]; then
    echo usage $0 filename
    exit 1
fi


filename="$1"


if [ ! -f $filename ]; then
    echo $filename does not exist
    exit 1
fi

echo $filename

set -x

paint_stars.jl distribution_function.hdf5 energies.hdf5 $filename

stellar_profiles_3d.jl ../../fiducial/out/ ${filename%".toml"}_stars.hdf5
