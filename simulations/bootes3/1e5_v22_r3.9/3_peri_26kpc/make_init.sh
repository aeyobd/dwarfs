#!/bin/bash

orbitpath=$DWARFS_ROOT/orbits/bootes3/EP2020_special_cases/orbit_peri_26
rm orbit.csv orbit.toml

set -xe
cp $orbitpath.csv orbit.csv
cp $orbitpath.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv > setup.log 2>&1
