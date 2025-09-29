
if [ $# -ne 1 ]; then
    echo "Usage: $0 starsname"
    exit 1
fi

set -xe 

out_path="../combined.hdf5"

iso_stars_path="../../stars/"
stars_file=$iso_stars_path/$1/probabilities_stars.hdf5
idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)

stellar_profile_3d.jl $out_path/$idx_f $stars_file $1/stellar_profile_f_3d.toml 
stellar_profile_3d.jl $out_path/$idx_f $stars_file $1/stellar_profile_unbound_f_3d.toml  -u
stellar_profile_3d.jl $out_path/1 $stars_file $1/stellar_profile_i_3d.toml 

