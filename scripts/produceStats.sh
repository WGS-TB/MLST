#!/bin/bash
sample=/home/glgan/Documents/Borrelia/Test_Data/SRR2034333_old/
files=`ls $sample`

for f in $files
do
	echo $f
	cd $sample$f
	python ~/Documents/Borrelia/scripts/Compute_values.py -g $f -l 41 > ~/Documents/Borrelia/$1/$f\_output_stats.txt -o $1
done	
