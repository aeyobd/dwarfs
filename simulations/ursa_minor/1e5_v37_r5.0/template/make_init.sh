#!/bin/bash

# rescales input relative to the halo
# TODO: change the halo to the correct model for the run
cp $DWARFS_ROOT/analysis/MODEL_PATH/ORBIT orbit.csv

julia $LGUYS_SCRIPTS/set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
