#!/bin/bash

orbitpath=$DWARFS_ROOT/analysis/sculptor/1e5_new_v43_r7/orbit_smallperi.1/next_orbit
rm orbit.csv orbit.toml
rm initial.hdf5

set -xe
cp $orbitpath.csv orbit.csv
cp $orbitpath.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv 
