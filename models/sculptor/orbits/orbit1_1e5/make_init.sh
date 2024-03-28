$LGUYS_SCRIPTS/centre_snapshot.jl -c isolation.hdf5 centred.hdf5
$LGUYS_SCRIPTS/set_in_orbit.jl centred.hdf5 initial.hdf5 -f orbit.csv
