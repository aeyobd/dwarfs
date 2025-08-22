#!/bin/bash

orbitpath=$DWARFS_ROOT/analysis/sculptor/1e6_new_v31_r3.2/L3M11_9Gyr_smallperi/orbit_v_R
rm orbit.csv orbit.toml
rm initial.hdf5

set -xe
cp $orbitpath.csv orbit.csv
cp $orbitpath.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
