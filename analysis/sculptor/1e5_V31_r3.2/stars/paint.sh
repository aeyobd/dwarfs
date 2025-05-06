#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage $0 dirname"
    exit 1
fi

paint_stars.jl distribution_function.hdf5 energies.hdf5 $1/profile.toml $1/probabilities --r-max 100 
