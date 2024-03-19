
$LGUYS_SCRIPTS/nbody.jl
$LGUYS_SCRIPTS/rescale.jl model.hdf5 isolation.hdf5 -m '20e-10' -r 0.0015
$LGUYS_SCRIPTS/centre_snapshot.jl isolation.hdf5 isolation.hdf5
$LGUYS_SCRIPTS/set_in_orbit.jl isolation.hdf5 initial.hdf5 -f orbit.csv
