#!/bin/bash

orbitpath=$DWARFS_ROOT/analysis/galaxy/simulation/orbit_name/orbit_base
rm orbit.csv orbit.toml

set -xe
cp $orbitpath.csv orbit.csv
cp $orbitpath.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv 
