#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo usage $0 galaxyname
    exit 1
fi

sbatch <<EOT
#!/bin/bash
#SBATCH -t 12:0:0
#SBATCH -n 16
#SBATCH -A durham
#SBATCH -p cosma5
#SBATCH -J $1.mcmc_2exp

#SBATCH -o "../"$1"/mcmc/mcmc_2exp.log"
export JULIA_NUM_THREADS=16
julia mcmc_2exp.jl $1

EOT

