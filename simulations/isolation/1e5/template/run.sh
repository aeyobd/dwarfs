#!/bin/bash 
mpirun -np $SLURM_NTASKS $GADGET_PATH/Gadget param.txt > log.out
