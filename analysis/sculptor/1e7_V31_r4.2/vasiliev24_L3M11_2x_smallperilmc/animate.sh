set -xe 

#project_potential.jl potential_stars.ini --limits 400 -n 1024 -o projected_mw_stars.hdf5 -T 0

# project density
#project_potential.jl lmc_potential.ini -T agama_times.txt -k 1 --limits 400 -n 1024 -o projected_lmc.hdf5
#project_potential.jl mw_halo_potential.ini -T agama_times.txt -k 1 --limits 400 -n 1024 -o projected_mw_halo.hdf5
#project_2d.jl . --limits 400 -k 1 -n 1024
#
t_f=$(tomlq .t_f_gyr orbital_properties.toml)

# make pictures (:)
files="projected_mw_halo.hdf5 projected_lmc.hdf5 projected_mw_stars.hdf5 projected_densities.hdf5 stars/exp2d_rs0.13/projected_densities.hdf5 "
colors="grey grey grey purple yellow"
scales="match2 0.5max match2 max max"

animate_dm.jl -i $files --scalebar 100 -s $scales --colors $colors -P 0.5 --time-today $t_f --time-scale 207.4

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

