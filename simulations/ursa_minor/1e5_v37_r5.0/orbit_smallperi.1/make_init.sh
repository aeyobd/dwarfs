#!/bin/bash

# rescales input relative to the halo
# TODO: change the halo to the correct model for the run
set -xe

old_analy_path=$DWARFS_ROOT/analysis/ursa_minor/1e5_v37_r5.0/orbit_smallperi
cp $old_analy_path/next_orbit.csv orbit.csv
cp $old_analy_path/next_orbit_shifts.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
