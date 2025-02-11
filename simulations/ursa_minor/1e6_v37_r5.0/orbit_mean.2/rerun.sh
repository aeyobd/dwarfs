#!/bin/bash 
mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMW param.2.txt 2 209 >> log.out
