#!/bin/bash

if [ "$#" -ne 1 ];
    then echo "usage: $0 rescale name"
    exit
fi


rescale_nfw.jl fiducial.hdf5 -o $1/isolation.hdf5 -n fiducial.toml -p $1/halo.toml 
