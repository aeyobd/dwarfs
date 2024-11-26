# make frames
set -xe 

animate_dm.jl -n 1001 -p ~/dwarfs/agama/potentials/vasiliev+21/potential_evolving.ini  --limits 600 --timescale 0.004822 -k 1 -P 0.5 -e

# combine to movie!!!
cd figures/dm_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4

