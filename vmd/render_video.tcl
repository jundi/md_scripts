# Number of frames
set num [molinfo top get numframes]

for {set i 0} {$i < $num} {incr i} {
	animate goto $i
	set filename snap_[format "%04d" $i].ppm
	render TachyonLOptiXInternal $filename
}
