#!/bin/bash
#SBATCH -A durham
#SBATCH -p cosma5
#SBATCH -t 6:0:0
#SBATCH --tasks 1
#SBATCH -n 1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu 2048

set -xe 

project_potential.jl ../potential_lmc_centre.ini -T times.txt -k 1 --limits 100 -n 257 -o projected_lmc.hdf5
