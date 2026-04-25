set -xe 

module load ffmpeg

cd animation
ffmpeg -framerate 50 -i frame_%d.png -c:v libx264 -crf 19 -pix_fmt yuv420p scl_animation.mp4

