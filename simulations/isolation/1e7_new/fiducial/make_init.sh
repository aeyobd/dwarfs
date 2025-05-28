#!/bin/bash
haloname="nfw_1e7_t20.0_xi3.0"
rescale_nfw.jl $DWARFS_ROOT/agama/halos/$haloname.hdf5 -n $DWARFS_ROOT/agama/halos/$haloname.toml -p halo.toml -o initial.hdf5 
