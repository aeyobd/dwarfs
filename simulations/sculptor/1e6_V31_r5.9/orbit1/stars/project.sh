#!/bin/bash

set -xe

#requires two arguments
if [ $# -ne 1 ]; then
    echo "Usage: ./project.sh <stars_file>"
    exit 1
fi

idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)

echo "using idx_f = $idx_f"


stars_path="../../../../isolation/1e6/fiducial/stars_ana/"
stars_file="$stars_path/$1_stars.hdf5"

if [ ! -f $stars_file ]; then
    echo "File $stars_file not found!"
    exit 1
fi

snap_path="../out/"

echo "writing stars"
project_snapshot.jl $snap_path $stars_file $1.fits -i $idx_f

echo "initial stars"
project_snapshot.jl $snap_path $stars_file $1_i.fits -i 1

# create time evolution
#
echo "making stellar profiles for each time"

stellar_profiles.jl .. $stars_file -o $1_stellar_profiles.hdf5
stellar_profiles_3d.jl .. $stars_file -o $1_stellar_profiles_3d.hdf5


