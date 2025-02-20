#!/bin/bash
haloname="nfw_1e4"
rescale_nfw.jl $DWARFS_ROOT/agama/halos/$haloname.hdf5 -n $DWARFS_ROOT/agama/halos/$haloname.toml -p halo.toml -o initial.hdf5 > rescale.log 2>&1
