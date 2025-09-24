#!/bin/bash
haloname="halos_cored/cnfw_1e5_rc_0.1"
rescale_nfw.jl $DWARFS_ROOT/agama/$haloname.hdf5 -n $DWARFS_ROOT/agama/$haloname.toml -p halo.toml -o initial.hdf5 
