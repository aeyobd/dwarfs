#!/bin/batch
#SBATCH --acount cosma
#SBATCH -p cosma5
#SBATCH --mem-per-cpu 4096
#SBATCH -c 1
#SBATCH -t 2:0:0
julia integrate.jl
