#!/bin/bash

source paths.sh

limits=150

proj_args="--x_vec 0 1 0 --y_vec 0 0 1 --limits $limits"


project_potential.jl $model_dir/potential_halo.ini $proj_args -T 0 -o projected_mw_halo.hdf5 
project_potential.jl $model_dir/potential_stars.ini $proj_args -T 0 -o projected_mw_stars.hdf5

project_2d.jl $model_dir $proj_args -k 1
project_2d.jl $model_dir $proj_args -k 1 -o projected_stars.hdf5 -s $stars_dir/probabilities_stars.hdf5
