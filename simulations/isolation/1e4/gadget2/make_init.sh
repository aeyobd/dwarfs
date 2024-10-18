#!/bin/bash
rescale_nfw.jl $DWARFS_ROOT/agama/halos/nfw_1e4.hdf5 -n $DWARFS_ROOT/agama/halos//halo.toml -p halo.toml -o initial_g4.hdf5
to_gadget2.jl initial_g4.hdf5 initial.hdf5
