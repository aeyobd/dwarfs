#!/bin/bash

# echo calculating energies
get_energies.jl ../out/20 energies.hdf5

echo calculating distribution function
calc_dist_func.jl energies.hdf5 distribution_function.hdf5 -n 100 -b equal_width 
