#!/bin/bash

rm orbit.csv
rm initial.hdf5
cp $DWARFS_ROOT/analysis/sculptor/mc_orbits/vasiliev24_L3M11_2x_special_cases/orbit_smallperilmc.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
