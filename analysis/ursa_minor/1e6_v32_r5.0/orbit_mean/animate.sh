set -x

# project density
rm projected_densities.hdf5 projected_mw.hdf5
rm figures/dm_animation/frame_*.png
rm figures/dm_animation/output.mp4

project_2d.jl . --limits 150 -k 1 -n 1001
project_potential.jl simulation/agama_potential.ini -T 0 --limits 150 -n 1001 -o projected_mw.hdf5

# make pictures (:)
#
agama_dir="$HOME/dwarfs/agama/potentials/animation/"
animate_dm.jl -i projected_densities.hdf5 projected_mw.hdf5 --scalebar 25 -s max max -P 0.5 -c yellow white 

cd figures/dm_animation
# combine to movie!!!
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

