
snap_path="$DWARFS_ROOT/analysis/isolation/1e7/fiducial/combined.hdf5/21"

halo_in="$DWARFS_ROOT/simulations/isolation/1e7/fiducial/halo.toml"
halo_out="../halo.toml"

set -ex


rescale_nfw.jl $snap_path -o initial_stars.hdf5 -n $halo_in -p $halo_out
get_energies.jl initial_stars.hdf5 energies.hdf5
calc_dist_ana.jl energies.hdf5 $halo_out distribution_function.hdf5
