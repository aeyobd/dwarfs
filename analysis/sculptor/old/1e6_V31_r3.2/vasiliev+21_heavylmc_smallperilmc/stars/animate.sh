set -xe 

# project density
# project_2d.jl .. -s ../../stars/$1/probabilities_stars.hdf5 -o $1/2d_densities.hdf5 --limits 500 -k 1

# make pictures (:)
animate_dm.jl -i $1/2d_densities.hdf5 ~/dwarfs/agama/potentials/vasiliev+21/animation/projected_disk.hdf5 ~/dwarfs/agama/potentials/vasiliev+21/animation/projected_lmc.hdf5 --scalebar 100 -s max 0.1max 0.01match2 -P -6 -o $1/animation/

# combine to movie!!!
cd $1/animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

