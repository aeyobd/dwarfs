set -xe 

# project density
#project_potential.jl simulation/agama_potential.ini -T agama_times.txt -k 1 --limits 500 -n 101 -o projected_lmc.hdf5
project_2d.jl . --limits 500 -k 1 -n 101

# make pictures (:)
animate_dm.jl -i projected_densities.hdf5 projected_lmc.hdf5 --scalebar 100 -s max max -P 0.5

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 20 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

