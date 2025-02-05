#!/bin/bash

rm orbit.csv
rm initial.hdf5
cp $DWARFS_ROOT/analysis/sculptor/mc_orbits/ORBIT_DIR/ORBIT.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
