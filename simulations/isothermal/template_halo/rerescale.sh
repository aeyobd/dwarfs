#!/bin/bash

source iso_paths.sh

rescale_nfw.jl $iso_output/$iso_idx_f -o isolation.hdf5 -n $iso_output/halo.toml -p halo.toml 
