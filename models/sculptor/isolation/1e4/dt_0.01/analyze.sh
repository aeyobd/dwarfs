#!/bin/bash

cd out
combine_outputs.py
calc_centres.jl --r_max 10
combine_outputs.py centres.hdf5
calc_profiles.jl -k 1
