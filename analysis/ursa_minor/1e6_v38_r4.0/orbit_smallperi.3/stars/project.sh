if [ $# -ne 1 ]; then
    echo "Usage: $0 starsname"
    exit 1
fi


mkdir -p $1

set -xe 

idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)
dist=$(awk -F' = ' '/distance_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)
source ../../paths.sh

iso_stars_path="../../stars/"
stars_file=$iso_stars_path/$1/probabilities_stars.hdf5

out_path="../combined.hdf5"
project_snapshot.jl $out_path $stars_file $1/final.fits -i $idx_f
project_snapshot.jl $out_path $stars_file $1/initial.fits -i 1 --distance $dist
