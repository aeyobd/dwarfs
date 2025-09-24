#!/bin/bash 
mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMW_AU param.txt > log.out
