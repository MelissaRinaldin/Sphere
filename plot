#!/bin/bash

file=`ls s-t*.dat | sort -n -k 1.4`

##############

zero=`ls s-t*.dat | head -n1` 
size=`wc -l $zero | awk '/dat/ {print $1}'`

num_of_lines=`echo "sqrt($size)" | bc`
box_size=`echo "0.078125*(${num_of_lines}-1)" | bc`

##############

numb=0

for f in $file; do

echo "Processing file $f ..."

name="${numb}.png"

if [ $numb -lt 1000 ]; then
	name="0${numb}.png"
fi

if [ $numb -lt 100 ]; then
	name="00${numb}.png"
fi

if [ $numb -lt 10 ]; then
	name="000${numb}.png"
fi

rm -f tmp*.dat

for (( n=1; n<=num_of_lines; n++ )); do
	
	head -n$((n*num_of_lines)) $f | tail -n$num_of_lines >> tmp.dat
	echo "" >> tmp.dat
	
done

command_line="'tmp.dat' us 1:2:3 notitle"

gnuplot <<EOF
set size square
set xrange [0:$box_size]
set yrange [0:$box_size]
set pm3d map
set palette gray
#set palette defined (-1 "red", 0 "white", 1 "green")
set term png crop medium size 800,800
unset tics
unset colorbox
set out '$name'
splot $command_line
EOF

numb=$((numb+1))

done

rm -f tmp*.dat

ffmpeg -r 7 -i %04d.png -sameq video.avi

rm -f *.png
