#!/bin/bash

#requires two arguments
if [ $# -ne 1 ]; then
    echo "Usage: ./project.sh <stars_file>"
    exit 1
fi

idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)

echo "using idx_f = $idx_f"


stars_path="../../../../isolation/1e6/fiducial/ana_stars/"
stars_file="$stars_path/$1_stars.hdf5"

if [ ! -f $stars_file ]; then
    echo "File $stars_file not found!"
    exit 1
fi

snap_path="../out/"

echo "writing stars"
project_snapshot.jl $snap_path $stars_file $1.fits -i $idx_f


# project initial stars too
idx_f=1
distance=81.4656 # take from orbital_properties.toml

echo "writing initial stars"
project_snapshot.jl $snap_path $stars_file $1_$idx_f.fits -i $idx_f -d $distance
