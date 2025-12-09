#!/bin/bash
stellar_profile.jl ../plummer_rs0.20/initial.fits final_profile.toml --bin-method both --centre-method first
stellar_profile.jl ../plummer_rs0.20/final.fits initial_profile.toml --bin-method both --centre-method first
