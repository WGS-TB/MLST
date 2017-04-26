#!/bin/bash

gene=$1 #gene name
sample=$2 #sample
reference=$3
#bowtie-build ./$gene".fas" $gene"_bowtie"		#build index of reference gene
#echo "Building index for "$gene" done."
#echo ""
#bowtie -a -p 8 $gene"_bowtie" -1 ./$sample"_1.fa" -2 ./$sample"_2.fa" $sample"_"$gene".out"	#mapping sample sequence to reference gene
#awk '{print $1"\t"$3"\t"$8}' $sample"_"$gene".out" > $sample"_"$gene"_reads.txt" #create table for python script
#echo "Mapping "$sample" sequence to "$gene" done."
#python /home/glgan/Documents/Borrelia/scripts/preprocess.py --path $sample'_'$gene'_reads.txt' --gene $gene #run the python script
#echo "computing proportions of "$gene" done"
bowtie -a -p 8 $reference -1 ./$sample"_1.fa" -2 ./$sample"_2.fa" $gene".out"  >/dev/null 2>&1
awk '{print $1"\t"$3"\t"$8"\t"$5}' $gene".out" > $gene"_reads.txt" 
rm ClpX_*.fas
echo ""

