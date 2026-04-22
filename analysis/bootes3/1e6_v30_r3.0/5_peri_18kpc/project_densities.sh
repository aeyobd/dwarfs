#!/bin/bash

limits=150

project_potential.jl potential_halo.ini --limits $limits -T 0 -o projected_mw_halo.hdf5
project_potential.jl potential_stars.ini --limits $limits -T 0 -o projected_mw_stars.hdf5

project_2d.jl . --limits $limits -k 1
project_2d.jl . --limits $limits -k 1 -o projected_stars.hdf5 -s ../stars/exp2d_rs0.20/probabilities_stars.hdf5
