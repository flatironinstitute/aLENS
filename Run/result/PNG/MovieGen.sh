#! /bin/bash
# 24 is the frame rate. recomment between 15 and 30
# crf controls quality. 17 is high, 24 is low, etc
# 1280 is width pixel count, -2 auto scales the height accordingly
# change moviename to the actual name

#get grandparent folder name

GRD="$(cd ../../; pwd)"
NAME1=$(basename $GRD)
GRD="$(pwd)"
NAME2=$(basename $GRD)

MOVIENAME=${NAME1}_${NAME2}
echo $MOVIENAME

ffmpeg -framerate 30 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p -profile:v high -level 4.2 -crf 17 -vf scale=1920:-2 -r 25 $MOVIENAME.mp4 
