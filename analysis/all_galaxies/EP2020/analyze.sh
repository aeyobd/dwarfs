#!/bin/bash

rm simulation
rm *.hdf5 
rm *.csv


simdir="$DWARFS_ROOT/simulations/all_galaxies/EP2020"
ln -s $simdir simulation

combine_outputs.py simulation/out
