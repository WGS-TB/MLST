#!/bin/bash
#This script takes 4 arguments: 
#1) path to python script i.e. borreliaILP.py
#2) path to directory containing data of each samples(can only contain reference.csv and directories of samples)
#3) path to reference MLST
#4) output path for conclusion
py_script_path=$1
data_path=$2
ref_path=$3
output_path=$4

mkdir $output_path
python $py_script_path --data $data_path --ref $ref_path --output $output_path 
