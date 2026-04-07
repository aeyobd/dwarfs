#!/bin/bash

if [ "$#" -ne 0 ]; then
    echo usage $0 
    exit 1
fi

sbatch <<EOT
#!/bin/bash
#SBATCH -t 3:0:0
#SBATCH -n 16
#SBATCH -A durham
#SBATCH -p cosma5
#SBATCH -J mcmc_gaia

#SBATCH -o "./mcmc/mcmc_gaia.log"
export JULIA_NUM_THREADS=16
julia mcmc_gaia.jl 

EOT

