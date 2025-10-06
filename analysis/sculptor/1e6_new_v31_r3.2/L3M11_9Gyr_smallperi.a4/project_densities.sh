#!/bin/bash
#SBATCH -p cosma5
#SBATCH -A durham
#SBATCH -t 2:0:0
#SBATCH --tasks 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 400

limits=400
bins=1024

#project_2d.jl . --limits $limits -k 1 -n $bins
project_potential.jl potential_stars.ini -o projected_stars.hdf5 -T 0 --limits $limits -n $bins

