
cp ../../../../isolation/1e6/halos/V60_r5.4/halo.toml .
cp ../../../../isolation/1e6/halos/V60_r5.4/isolation.hdf5 .

$LGUYS_SCRIPTS/set_in_orbit.jl isolation.hdf5 initial.hdf5 -f ../orbit.csv
