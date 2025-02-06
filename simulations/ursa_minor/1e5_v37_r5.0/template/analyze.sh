#!/bin/bash

set -e 

cd out
combine_outputs.py
calc_centres.jl --r_max 1
combine_outputs.py centres.hdf5
mass_profiles.jl 
