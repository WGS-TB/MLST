#!/bin/bash
sample=$1
sample_path=$2
output_path=$3
var_db=$4

loci="clpA clpX nifS pepX pyrG recG rplB uvrA"

cd $output_path
mkdir $sample
cd $sample
echo "now processing "$sample
for l_i in $loci
do
	bowtie-build $var_db"/"$l_i".fas" $l_i"_bowtie" >/dev/null
	echo ""
	echo "building index for "$l_i" done"
	echo ""
	#map the sample to the locus
	echo "Mapping sample to variant....."
	bowtie -a -p 8 $l_i"_bowtie" -1 $sample_path"/"$sample"_1.fastq" -2 $sample_path"/"$sample"_2.fastq" $sample"_"$l_i".out" >/dev/null 
	#create reads table
	echo "Creating reads table....."
	awk '{print $1"\t"$5"\t"$10}' $sample"_"$l_i".out" > $sample"_"$l_i"_reads.txt"
	#generate matrix, predict variants, and compute proportions
	echo "predicting variants and computing proportions........"
	python ~/Desktop/SRA_bowtie/Alignments/scripts/pipeline_wrapper.py -i $sample"_"$l_i"_reads.txt" -g $l_i
done
rm *.ebwt
rm *.out