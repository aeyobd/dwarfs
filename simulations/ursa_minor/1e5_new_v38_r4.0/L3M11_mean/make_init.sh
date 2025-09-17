#!/bin/bash

orbitpath=$DWARFS_ROOT/orbits/ursa_minor/vasiliev24_L3M11_special_cases/orbit_mean
rm orbit.csv orbit.toml

set -xe
cp $orbitpath.csv orbit.csv
cp $orbitpath.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv 
