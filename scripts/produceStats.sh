#!/bin/bash
sample=/home/glgan/Documents/Borrelia/Test_Data/SRR2034333_old/
files=`ls $sample`
outputFolderPath=$1

for f in $files
do
	echo $f
	cd $sample$f
	rm *.out *reads.txt *.fa
	python ~/Documents/Borrelia/scripts/Compute_values.py -g $f -l 40 -o $outputFolderPath > $outputFolderPath$f\_output_stats.txt 
done	
