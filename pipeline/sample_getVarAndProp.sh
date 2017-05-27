#!/bin/bash
#This script takes 4 arguments: Sample name, path to sample reads, output path for identified variants and their proportions, indexed reference variants for each loci
#Finally it will run getVariantAndProp.py which outputs identified variants and their proportions
sample=$1
sample_path=$2
output_path=$3
var_db=$4
py_script_path=$5

loci="clpA clpX nifS pepX pyrG recG rplB uvrA"

cd $output_path

if [ ! -d $sample ]
then
	mkdir $sample
fi

cd $sample
echo ""
echo "============================================================================================="
echo "==============================   Now processing "$sample"   ================================"
echo "============================================================================================="
for l_i in $loci
do
	echo ""
	echo "~~~~~~~~~~~~~~~~~~~~~~~ Gene "$l_i" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	echo ""
	#building index for gene
	bowtie-build $var_db"/"$l_i".fas" $l_i"_bowtie" >/dev/null
	echo ".... Done building index ...."
	#map the sample to the locus
	echo ".... Mapping sample reads to variants ...."
	echo ""
	bowtie -a -v 3 -p 8 $l_i"_bowtie" -1 $sample_path"/"$sample"_1.fastq" -2 $sample_path"/"$sample"_2.fastq" $sample"_"$l_i".out" >/dev/null 
	#create reads table
	echo ".... Summarizing Bowtie mapping information into reads table ....."
	echo ""
	awk '{print $1"\t"$5"\t"$10}' $sample"_"$l_i".out" > $sample"_"$l_i"_reads.txt"
	#generate matrix, predict variants, and compute proportions
	echo "....... Predicting variants and computing proportions ........"
	echo ""
	python $py_script_path -i $sample"_"$l_i"_reads.txt" -g $l_i
	echo ""
done
rm *.ebwt
rm *.out
