#!/bin/bash
 
rm -f out.mp4
cd fluid
convert '*.png[640x]' resized_%03d.jpg
cd ..
ffmpeg -i fluid/resized_%3d.jpg -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
