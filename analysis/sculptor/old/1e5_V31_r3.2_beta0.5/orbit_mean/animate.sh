set -xe 

# project density
project_2d.jl . --limits 200 -k 1

# make pictures (:)
#
agama_dir="$HOME/dwarfs/agama/potentials/animation/"
animate_dm.jl -i projected_densities.hdf5 $agama_dir/projected_all.hdf5 --scalebar 100 -s max max -P 0.5 -c yellow white 

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

