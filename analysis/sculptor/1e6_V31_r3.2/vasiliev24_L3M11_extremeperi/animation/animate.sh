set -xe 

# make pictures (:)
animate_dm.jl -i projected_densities.hdf5 projected_mw.hdf5 projected_lmc.hdf5 --scalebar 25 -s max max max -P 0.5 -c yellow white orange

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 20 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

