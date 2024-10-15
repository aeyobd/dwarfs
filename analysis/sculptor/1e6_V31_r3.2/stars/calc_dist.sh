source ../paths.sh

halo_out="../halo.toml"

set -ex
snap_path="$iso_out/$iso_idx_f"


rescale_nfw.jl $snap_path -o initial_stars.hdf5 -n $halo_in -p $halo_out
get_energies.jl initial_stars.hdf5 energies.hdf5
calc_dist_ana.jl energies.hdf5 $halo_out distribution_function.hdf5
