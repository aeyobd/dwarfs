#!/bin/bash 
#restartoption="1"
#restartoption="2 snapnum"
mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMW param.txt $restartoption >> log.out
