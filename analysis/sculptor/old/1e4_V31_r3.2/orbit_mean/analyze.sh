#!/bin/bash

source paths.sh

set -xe

cp $out_path/agama_potential.ini .
exit 0
combine_outputs.py $out_path/out/
calc_centres.jl --r_max 10
combine_outputs.py $out_path/out -c centres.hdf5
mass_profiles.jl -k 1
