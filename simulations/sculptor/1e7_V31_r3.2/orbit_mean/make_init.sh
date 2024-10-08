#!/bin/bash

cp $DWARFS_ROOT/analysis//sculptor/mc_orbits/orbit1.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
