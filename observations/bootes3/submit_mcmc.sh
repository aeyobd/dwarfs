#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo usage $0 samplename
    exit 1
fi

sbatch <<EOT
#!/bin/bash
#SBATCH -t 5:0:0
#SBATCH -n 16
#SBATCH -A durham
#SBATCH -p cosma5
#SBATCH -J mcmc/$1.mcmc_plummer

#SBATCH -o "./mcmc/$1.mcmc_plummer.log"
export JULIA_NUM_THREADS=16
julia mcmc_mf_prof.jl $1 

EOT

