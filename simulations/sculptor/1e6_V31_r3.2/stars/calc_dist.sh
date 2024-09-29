#!/bin/bash

echo calculating energies
get_energies.jl ../mid.hdf5 energies.hdf5

echo calculating distribution function
calc_dist_ana.jl energies.hdf5 ../halo.toml distribution_function.hdf5
