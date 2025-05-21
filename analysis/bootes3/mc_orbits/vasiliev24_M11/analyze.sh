#!/bin/bash
out_dir=$DWARFS_ROOT/simulations/bootes3/mc_orbits/vasiliev24_M11/
rm simulation
rm *.hdf5 *.fits

ln -s $out_dir simulation

combine_outputs.py $out_dir/out
julia ../calc_peris.jl
