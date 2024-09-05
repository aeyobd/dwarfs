# rescales input relative to the halo
julia $LGUYS_SCRIPTS/rescale_nfw.jl $LGUYS_SCRIPTS/../agama/halos/nfw_1e6.hdf5 -n $LGUYS_SCRIPTS/../agama/halos/halo.toml -p halo.toml -o initial.hdf5
