#!/bin/bash

set -xe 

cp $DWARFS_ROOT/analysis/sculptor/mc_orbits/vasiliev_lmc_special_cases/orbit_smallperi.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
