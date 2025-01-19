#!/bin/bash
rescale_nfw.jl $DWARFS_ROOT/agama/halos/asy_1e6.hdf5 -n $DWARFS_ROOT/agama/halos/halo_asy.toml -p halo.toml -o initial.hdf5
