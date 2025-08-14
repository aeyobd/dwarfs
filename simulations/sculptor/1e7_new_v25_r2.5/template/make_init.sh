#!/bin/bash

potentialpath=$DWARFS_ROOT/orbits/sculptor/potentialname
orbitname=orbit_
rm orbit.csv orbit.toml
rm initial.hdf5

set -xe
cp $potentialpath/$orbitpath.csv orbit.csv
cp $potentialpath/$orbitpath.toml orbit.toml
cp $potentialpath/agama_potential.ini .

set_in_orbit.jl ../isolation.hdf5 initial.hdf5 -f orbit.csv
