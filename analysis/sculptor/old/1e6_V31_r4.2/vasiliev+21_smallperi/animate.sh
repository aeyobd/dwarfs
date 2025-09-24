set -xe 

# project density
# project_2d.jl . --limits 500 -k 1

# make pictures (:)
#
agama_dir="$HOME/dwarfs/agama/potentials/vasiliev+21/animation"
animate_dm.jl -i projected_densities.hdf5 $agama_dir/projected_halo.hdf5 $agama_dir/projected_lmc.hdf5 $agama_dir/projected_disk.hdf5 --scalebar 100 -s max match4 match4 max -P 0.5 -c yellow white orange white

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

