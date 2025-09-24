set -xe 

# project density

project_2d.jl ../.. --limits 400 -n 1024 -k 1 -s ../../../stars/exp2d_rs0.13/probabilities_stars.hdf5

