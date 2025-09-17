#!/bin/bash
#SBATCH --account durham
#SBATCH --partition cosma5
#SBATCH --mem-per-cpu 4096
#SBATCH --cpus-per-task 1
#SBATCH --time 4:0:0
#SBATCH --job-name "L3M11 mcmc orbits"

rm -rf out/*
julia integrate.jl
