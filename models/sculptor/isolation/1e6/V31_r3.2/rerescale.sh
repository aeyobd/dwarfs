#!/bin/bash

rescale_nfw.jl ../fiducial/mid.hdf5 -o mid.hdf5 -n ../fiducial/halo.toml -p halo.toml 
rescale_nfw.jl ../fiducial/final.hdf5 -o final.hdf5 -n ../fiducial/halo.toml -p halo.toml 
