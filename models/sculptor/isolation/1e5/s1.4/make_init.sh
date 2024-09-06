#!/bin/bash

# rescales input relative to the halo
# TODO: change the halo to the correct model for the run
julia $LGUYS_SCRIPTS/rescale_nfw.jl $LGUYS_SCRIPTS/../agama/halos/nfw_1e5.hdf5 -n $LGUYS_SCRIPTS/../zeno/halo.toml -p halo.toml -o initial.hdf5
