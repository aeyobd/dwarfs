set -xe 

#project_potential.jl potential_stars.ini --limits 400 -n 1024 -o projected_mw_stars.hdf5 -T 0

# project density
#project_potential.jl lmc_potential.ini -T agama_times.txt -k 1 --limits 400 -n 1024 -o projected_lmc.hdf5
#project_potential.jl mw_halo_potential.ini -T agama_times.txt -k 1 --limits 400 -n 1024 -o projected_mw_halo.hdf5
#project_2d.jl . --limits 400 -k 1 -n 1024
#
t_f=$(tomlq .t_f_gyr orbital_properties.toml)

# make pictures (:)
files="projected_lmc_resampled.hdf5"

animate_dm.jl -i $files --scalebar 100 -P 0.5 --time-today $t_f

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

