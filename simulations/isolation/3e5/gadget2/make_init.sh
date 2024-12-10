#!/bin/bash
ln -s ../fiducial/initial.hdf5 initial_g4.hdf5
to_gadget2.jl initial_g4.hdf5 initial.hdf5
