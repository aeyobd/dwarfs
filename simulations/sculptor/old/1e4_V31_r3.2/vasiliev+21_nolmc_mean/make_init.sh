#!/bin/bash

cp $DWARFS_ROOT/analysis/sculptor/mc_orbits/vasiliev_nolmc_special_cases/orbit_mean.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
