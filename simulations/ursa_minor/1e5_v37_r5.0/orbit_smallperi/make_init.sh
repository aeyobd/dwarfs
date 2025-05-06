#!/bin/bash

# rescales input relative to the halo
# TODO: change the halo to the correct model for the run
set -xe
export orbit_path=$DWARFS_ROOT/analysis/ursa_minor/mc_orbits/EP20_special_cases/orbit_smallperi

cp $orbit_path.csv orbit.csv
cp $orbit_path.toml orbit.toml

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
