set -xe 

# project density
# project_2d.jl . --limits 500 -k 1

# make pictures (:)
animate_dm.jl -i projected_densities.hdf5 --scalebar 100 -s max -P 0.5

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
