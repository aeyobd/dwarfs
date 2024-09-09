#!/bin/bash

#requires two arguments
if [ $# -lt 1 ]; then
    echo "Usage: ./project.sh stars_file [idx_f]"
    exit 0
fi

if [ $# -eq 2 ]; then
    idx_f=$2
else
    idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)
fi

echo "using idx_f = $idx_f"


stars_path="../../../../../isolation/1e6/halos/V60_r5.4/stars/"
stars_file="$stars_path/$1.hdf5"

if [ ! -f $stars_file ]; then
    echo "File $stars_file not found!"
    exit 1
fi

snap_path="../out/"

echo "writing stars"
project_snapshot.jl $snap_path $stars_file $1.fits -i $idx_f
