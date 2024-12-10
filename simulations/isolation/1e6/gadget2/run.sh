#!/bin/bash 
time  mpirun -np $SLURM_NTASKS $HOME/source/Gadget-2.0.7/Gadget2/Gadget2 param.txt > log.out
