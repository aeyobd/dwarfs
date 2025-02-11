#!/bin/bash
rescale_nfw.jl $DWARFS_ROOT/agama/halos/nfw_1e5_beta0.5.hdf5 -n $DWARFS_ROOT/agama/halos//halo.toml -p halo.toml -o initial.hdf5
