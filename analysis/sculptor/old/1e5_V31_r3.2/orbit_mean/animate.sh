set -xe 

# project density
project_potential.jl simulation/agama_potential.ini --limits 200 -T 0 -o projected_potential.hdf5
project_2d.jl . --limits 200 -k 1

# make pictures (:)
#
animate_dm.jl -i projected_densities.hdf5 projected_potential.hdf5 --scalebar 100 -s max max -P 0.5 -c yellow white 

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 10 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

