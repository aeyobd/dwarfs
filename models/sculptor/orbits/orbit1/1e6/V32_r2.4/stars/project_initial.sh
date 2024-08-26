#!/bin/bash

#requires two arguments
if [ $# -ne 2 ]; then
    echo "Usage: ./project.sh stars_file [idx_f]"
    exit 0
fi

idx_f=$2

distance=77.8284
echo "using idx_f = $idx_f"


stars_path="../../../../../isolation/1e6/halos/V32_r2.4/stars/"
stars_file="$stars_path/$1_stars.hdf5"

if [ ! -f $stars_file ]; then
    echo "File $stars_file not found!"
    exit 1
fi

snap_path="../out/"

echo "writing stars"
project_snapshot.jl $snap_path $stars_file $1_$2.fits -i $idx_f -d $distance
