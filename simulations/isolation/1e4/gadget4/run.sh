#!/bin/bash 
time mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMono param.txt > log.out
