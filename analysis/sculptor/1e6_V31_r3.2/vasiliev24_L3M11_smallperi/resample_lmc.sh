source paths.sh
lmc_file="$out_path/trajlmc.txt"
resample_orbit.jl $lmc_file -o lmc_traj.csv -t combined.hdf5 
resample_orbit.jl $lmc_file -o lmc_traj_exp.csv -t orbit.csv
