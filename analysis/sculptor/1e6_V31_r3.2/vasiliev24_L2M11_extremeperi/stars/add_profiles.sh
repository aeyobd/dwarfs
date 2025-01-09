
if [ $# -ne 1 ]; then
    echo "Usage: $0 starsname"
    exit 1
fi


set -xe 

out_path="../combined.hdf5"

iso_stars_path="../../stars/"
stars_file=$iso_stars_path/$1/probabilities_stars.hdf5

stellar_profiles_3d.jl $out_path $stars_file -o $1/stellar_profiles_3d.hdf5 
#stellar_profiles.jl $out_path $stars_file -o $1/stellar_profiles.hdf5 --bin-method both

