#!/bin/bash
haloname="../halos_elliptical/1e6_oblate_0.5"
rescale_nfw.jl $DWARFS_ROOT/agama/halos/$haloname.hdf5 -n $DWARFS_ROOT/agama/halos/halo.toml -p halo.toml -o initial.hdf5 
