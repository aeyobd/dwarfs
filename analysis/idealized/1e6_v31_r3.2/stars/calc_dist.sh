source ../paths.sh

halo_out="../halo.toml"

set -ex

rescale_nfw.jl $iso_path/$iso_idx_paint -o iso_paint.hdf5 -n $iso_halo -p $halo_out > rescale.log 2>&1
rescale_nfw.jl $iso_path/1 -o iso_initial.hdf5 -n $iso_halo -p $halo_out >> rescale.log 2>&1
rescale_nfw.jl $iso_path/$iso_idx_f -o iso_final.hdf5 -n $iso_halo -p $halo_out >> rescale.log 2>&1

get_energies.jl iso_paint.hdf5 energies.hdf5 > energies.log 2>&1
calc_dist_ana.jl energies.hdf5 $halo_out distribution_function.hdf5 > dist.log 2>&1
