#!/bin/bash

sample=/home/glgan/Documents/Borrelia/Test_Data/SRR2034333/
files=`ls $sample`

for f in $files
do
	cd $sample$f
	python ~/Documents/Borrelia/scripts/Compute_values.py -g $f -l 41 > ~/Documents/Borrelia/simulated_stats/$f\_output_stats.txt
done	
