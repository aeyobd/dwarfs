#!/bin/bash
rescale_nfw.jl $DWARFS_ROOT/agama/halos/asy_1e4.hdf5 -n $DWARFS_ROOT/agama/halos/halo_asy.toml -p halo.toml -o initial.hdf5
