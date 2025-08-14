#!/bin/bash

calc_centres.jl --r_max 0.3 -o centres_r0.3.hdf5
calc_centres.jl -m COM -o centres_com.hdf5
calc_centres.jl -m MostBound -o centres_mb.hdf5
calc_centres.jl -m Potential -o centres_potential.hdf5
