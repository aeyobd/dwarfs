set -xe 


t_f=$(tomlq .t_f_gyr orbital_properties.toml)
idx_f=$(tomlq .idx_f orbital_properties.toml)

files="projected_mw_stars.hdf5 projected_stars.hdf5"
scales="0.001max max"
colors="grey yellow"

animate_dm.jl -i $files --scalebar 50 -s $scales -P -5 -c $colors --time-today $t_f --idx-max $idx_f -o figures/stars_animation

# combine to movie!!!
cd figures/stars_animation
ffmpeg -framerate 20 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

