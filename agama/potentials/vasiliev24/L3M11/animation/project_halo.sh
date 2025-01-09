
set -xe 

project_potential.jl ../potential_mw_halo.ini -T times.txt -k 1 --limits 500 -n 1001 -o projected_halo.hdf5
