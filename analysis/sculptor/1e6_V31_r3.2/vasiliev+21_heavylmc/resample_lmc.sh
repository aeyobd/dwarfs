lmc_file="$DWARFS_ROOT/simulations/sculptor/1e6_V31_r3.2/vasiliev+21_heavylmc_smallperilmc/trajlmc.txt"
resample_orbit.jl $lmc_file -o lmc_traj.csv -t combined.hdf5
resample_orbit.jl $lmc_file -o lmc_traj_exp.csv -t orbit.csv
