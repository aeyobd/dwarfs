#!/bin/bash

source paths.sh

set -xe

combine_outputs.py $out_path/out/
calc_centres.jl --r_max 0.3
combine_outputs.py $out_path/out -c centres.hdf5
mass_profiles.jl -k 1
