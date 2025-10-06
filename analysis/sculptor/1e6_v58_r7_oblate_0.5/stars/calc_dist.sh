source ../paths.sh

halo_out="../halo.toml"

set -ex

rescale_nfw.jl $iso_path/$iso_idx_paint -o iso_paint.hdf5 -n $iso_halo -p $halo_out
rescale_nfw.jl $iso_path/1 -o iso_initial.hdf5 -n $iso_halo -p $halo_out
rescale_nfw.jl $iso_path/$iso_idx_f -o iso_final.hdf5 -n $iso_halo -p $halo_out

#get_energies.jl iso_paint.hdf5 energies.hdf5
#calc_dist_func.jl energies.hdf5 distribution_function.hdf5
