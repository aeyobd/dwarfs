#!/bin/bash

# rescales input relative to the halo
# TODO: change the halo to the correct model for the run
cp ~/sculptor/isolation/1e6/V31_r3.2/final.hdf5 isolation.hdf5
cp ~/sculptor/isolation/1e6/V31_r3.2/halo.toml .
cp ~/sculptor/mc_orbits/orbit1.csv orbit.csv

julia $LGUYS_SCRIPTS/set_in_orbit.jl isolation.hdf5 initial.hdf5 -f orbit.csv
