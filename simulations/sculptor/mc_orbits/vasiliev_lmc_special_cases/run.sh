
#!/bin/bash 
mpirun -np $SLURM_NTASKS $GADGET_PATH/GadgetMC_AU param.txt > log.out
