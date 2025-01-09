#!/bin/bash

# remove analysis files
rm -rf figures
rm peris_apos.fits
rm combined.hdf5
rm orbit*.csv
rm orbit*.toml
rm trajectory.hdf5
rm lmc_traj.csv
rm initial.hdf5


out_dir=$DWARFS_ROOT/simulations/sculptor/mc_orbits/vasiliev_lmc_loose/out
combine_outputs.py $out_dir
julia resample_lmc.jl
julia calc_peris.jl
cp $out_dir/../initial.hdf5 .

mkdir figures
