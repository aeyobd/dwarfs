#!/bin/bash

# rescales input relative to the halo
# TODO: change the halo to the correct model for the run
cp $DWARFS_ROOT/analysis/ursa_minor/1e6_v37_r5.0/orbit_mean.1/new_orbit.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
