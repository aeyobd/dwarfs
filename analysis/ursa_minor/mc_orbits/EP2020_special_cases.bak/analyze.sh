#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/ursa_minor/mc_orbits/EP2020_special_cases/
rm simulation
rm *.hdf5

ln -s $out_dir simulation

combine_outputs.py simulation/out
julia ../calc_peris.jl
