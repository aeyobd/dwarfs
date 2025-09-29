#!/bin/bash
haloname="exp2d_1e4"
rescale_nfw.jl $DWARFS_ROOT/agama/halos_exp/$haloname.hdf5 -n $DWARFS_ROOT/agama/halos_exp/$haloname.toml -p halo.toml -o initial.hdf5 
