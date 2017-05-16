#!/bin/bash
sample=/home/glgan/Documents/Borrelia/Test_Data/SRR2034333_old/
files=`ls $sample`

for f in $files
do
	echo $f
	cd $sample$f
	python ~/Documents/Borrelia/scripts/Compute_values.py -g $f -l 40 -o $1 > ~/Documents/Borrelia/$1/$f\_output_stats.txt 
done	
