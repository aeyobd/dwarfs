
set -xe 

idx_i=103
dist=$(awk -F' = ' '/distance_f/ {gsub(/"/, "", $2); print $2}' ../../orbital_properties.toml)

source ../../../paths.sh

iso_stars_path="../../../stars/"

starsname=exp2d_rs0.10
stars_file=$iso_stars_path/$starsname/probabilities_stars.hdf5

out_path="../../combined.hdf5"
project_snapshot.jl $out_path $stars_file intermediate.fits -i $idx_i

