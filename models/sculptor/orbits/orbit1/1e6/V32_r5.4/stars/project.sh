#!/bin/bash

#requires two arguments
if [ $# -ne 1 ]; then
    echo "Usage: ./project.sh <stars_file>"
    exit 1
fi

idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)

echo "using idx_f = $idx_f"


stars_path="../../../../../isolation/1e6/stars/"
stars_file="$stars_path/$1.hdf5"

if [ ! -f $stars_file ]; then
    echo "File $stars_file not found!"
    exit 1
fi

snap_path="../out/"

echo "writing stars"
project_snapshot.jl $snap_path $stars_file $1.fits -i $idx_f