#!/bin/bash

orbitpath=$DWARFS_ROOT/analysis/ursa_minor/1e5_new_v38_r4.0/orbit_smallperi.2/next_orbit
rm orbit.csv orbit.toml

set -xe
cp $orbitpath.csv orbit.csv
cp $orbitpath.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv 
