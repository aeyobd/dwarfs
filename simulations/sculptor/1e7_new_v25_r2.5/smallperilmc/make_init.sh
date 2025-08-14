#!/bin/bash

potentialpath=$DWARFS_ROOT/orbits/sculptor/vasiliev24_L3M11_2x_special_cases
orbitname=orbit_smallperilmc
rm orbit.csv orbit.toml
rm initial.hdf5

set -xe
cp $potentialpath/$orbitname.csv orbit.csv
cp $potentialpath/$orbitname.toml orbit.toml
cp $potentialpath/agama_potential.ini .

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
