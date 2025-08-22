#!/bin/bash
#SBATCH -n 1
#SBATCH -A durham
#SBATCH -p cosma5
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4096
#SBATCH -t 16:0:0

export JULIA_NUM_THREADS=16
mass_profiles.jl -k 1
