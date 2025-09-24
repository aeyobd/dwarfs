set -xe 


t_f=$(tomlq .t_f_gyr orbital_properties.toml)
idx_f=$(tomlq .idx_f orbital_properties.toml)

files="projected_halo.hdf5 projected_stars.hdf5 projected_densities.hdf5"
scales="match2 0.5max max"
colors="grey grey purple"

animate_dm.jl -i $files --scalebar 50 -s $scales -P 0.5 -c $colors --time-today $t_f --idx-max $idx_f

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 20 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

