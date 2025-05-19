#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/vasiliev24_M11/
rm simulation
rm *.hdf5

ln -s $out_dir simulation

combine_outputs.py simulation/out
julia ../calc_peris.jl
