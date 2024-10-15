
if [ $# -ne 1 ]; then
    echo "Usage: $0 starsname"
    exit 1
fi


mkdir -p $1

set -xe 

idx_f=$(awk -F' = ' '/idx_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)
dist=$(awk -F' = ' '/distance_f/ {gsub(/"/, "", $2); print $2}' ../orbital_properties.toml)

out_path="../combined.hdf5"

iso_stars_path="../../stars/"
stars_file=$iso_stars_path/$1/probabilities_stars.hdf5


#project_snapshot.jl $out_path $stars_file $1/final.fits -i $idx_f
#project_snapshot.jl $out_path $stars_file $1/initial.fits -i 1 --distance $dist


source ../../paths.sh
out_path=$iso_path
project_snapshot.jl $out_path $stars_file --scale ../../halo-used.toml -i 1 --distance $dist $1/iso_initial.fits

project_snapshot.jl $out_path $stars_file --scale ../../halo-used.toml -i $idx_paint --distance $dist $1/iso_paint.fits

project_snapshot.jl $out_path $stars_file --scale ../../halo-used.toml -i $idx_f_iso --distance $dist $1/iso_final.fits


for f in $1/iso_initial.fits $1/iso_paint.fits $1/iso_final.fits; do
    stellar_profile.jl $f --mass-column probability -s --bin-method both
done
