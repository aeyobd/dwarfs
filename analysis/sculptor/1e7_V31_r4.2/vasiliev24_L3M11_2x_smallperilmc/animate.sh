set -xe 

# project density
#project_potential.jl lmc_potential.ini -T agama_times.txt -k 1 --limits 500 -n 1025 -o projected_lmc.hdf5
#project_potential.jl mw_halo_potential.ini -T agama_times.txt -k 1 --limits 500 -n 1025 -o projected_mw_halo.hdf5
project_2d.jl . --limits 500 -k 1 -n 1025

# make pictures (:)
animate_dm.jl -i projected_lmc.hdf5 projected_mw_halo.hdf5 projected.hdf5 
--scalebar 100 -s match2 max max -P 0.5

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 20 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

