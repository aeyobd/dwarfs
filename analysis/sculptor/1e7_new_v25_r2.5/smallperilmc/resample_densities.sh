
rm projected_lmc_resampled.hdf5 projected_mw_halo_resampled.hdf5

resample_density.jl projected_mw_halo.hdf5 projected_densities.hdf5 -o projected_mw_halo_resampled.hdf5 --time-scale 207.4


resample_density.jl projected_lmc.hdf5 projected_densities.hdf5 -o projected_lmc_resampled.hdf5 --time-scale 207.4 --centre-file trajlmc.txt
