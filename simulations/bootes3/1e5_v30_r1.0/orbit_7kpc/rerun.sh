#!/bin/bash 
mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMW param.txt 2 180 >> log.out
