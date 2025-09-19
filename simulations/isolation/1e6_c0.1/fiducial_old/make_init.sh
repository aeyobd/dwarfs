#!/bin/bash
rescale_nfw.jl $DWARFS_ROOT/agama/halos_cored/cnfw_1e6_rc_0.1.hdf5 -n $DWARFS_ROOT/agama/halos_cored/cnfw_1e6_rc_0.1.toml -p halo.toml -o initial.hdf5
