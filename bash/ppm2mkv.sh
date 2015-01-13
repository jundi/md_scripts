#!/bin/bash
#
# Creates high quality video from ppm images.
#

tune=animation		# film,animation,grain,stillimage,psnr,ssim
preset=slow		# ultrafast,superfast,veryfast,faster,fast,medium,slow,slower,veryslow,placebo
crf=20			# quality: 0-51 (smaller is better) 
framerate=50	
videoframerate=50
output=out.mkv

ffmpeg -framerate $framerate -pattern_type glob -i '*.ppm' -c:v libx264 -r $videoframerate -tune $tune -preset $preset -crf $crf $output
