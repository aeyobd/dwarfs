#!/bin/bash

# rescales input relative to the halo
cp $DWARFS_ROOT/analysis/sculptor/mc_orbits/ORBITMODEL/ORBIT_CAS.csv orbit.csv

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
