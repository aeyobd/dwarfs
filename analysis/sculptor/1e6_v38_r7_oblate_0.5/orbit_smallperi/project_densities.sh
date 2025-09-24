#!/bin/bash

limits=125

project_potential.jl potential_halo.ini --limits $limits -T 0 -o projected_halo.hdf5
project_potential.jl potential_stars.ini --limits $limits -T 0 -o projected_stars.hdf5

#project_2d.jl . --limits $limits -k 1
