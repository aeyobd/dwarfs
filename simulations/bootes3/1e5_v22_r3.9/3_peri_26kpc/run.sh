#!/bin/bash 
mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMW param.txt > log.out
