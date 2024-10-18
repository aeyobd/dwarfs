#!/bin/bash 
time  mpirun -np $SLURM_NTASKS $DWARFS_ROOT/gadget2_iso/Gadget2 param.txt > log.out
