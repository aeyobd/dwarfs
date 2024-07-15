# rescales input relative to fornax
julia $LGUYS_SCRIPTS/rescale_nfw.jl $LGUYS_SCRIPTS/../zeno/models/nfw_1e5.hdf5 -o initial.hdf5 -p halo.toml -n $LGUYS_SCRIPTS/../zeno/halo.toml
