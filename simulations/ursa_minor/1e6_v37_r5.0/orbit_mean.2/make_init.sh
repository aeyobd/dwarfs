#!/bin/bash


set -xe
cp $DWARFS_ROOT/analysis/ursa_minor/1e5_v37_r5.0/orbit_mean.2/orbit.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
