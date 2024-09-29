#!/bin/bash

nbody.jl exponential.hdf5
rescale.jl exponential.hdf5 initial.hdf5 -m 5e-6 -r 0.0107
