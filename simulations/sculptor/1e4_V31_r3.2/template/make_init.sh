#!/bin/bash

cp $DWARFS_ROOT/analysis/mc_orbits/sculptor/orbit_here orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
