
echo need to adjust which input halo to use!
exit 1

isolation_dir=../../../../../isolation/1e6
halo_name= # FILL IN!

cp $isolation_dir/halos/$halo_name/halo.toml .
cp $isolation_dir/halos/$halo_name/isolation.hdf5 .

$LGUYS_SCRIPTS/set_in_orbit.jl isolation.hdf5 initial.hdf5 -f ../orbit.csv
