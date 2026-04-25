set -xe 

module load ffmpeg

cd animation
ffmpeg -framerate 60 -i frame_%d.png -c:v libx264 -crf 19 -pix_fmt yuv420p scl_lmc_animation.mp4

