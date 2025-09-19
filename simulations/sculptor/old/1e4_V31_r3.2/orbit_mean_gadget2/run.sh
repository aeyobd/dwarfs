#!/bin/bash 
mpirun -np $SLURM_NTASKS $DWARFS_ROOT/gadget2_rapha/Gadget2-UVIC param.txt > log.out
