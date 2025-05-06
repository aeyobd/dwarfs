
#!/bin/bash 
mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMC param.txt > log.out
