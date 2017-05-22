#!/bin/bash
#1) data_path = path to directory which contains sample folders, where this folder contains reads info
#2) lociDB_path = path to loci databases folder
#3) ref_strains = path to reference strains
data_path=$1
lociDb_path=$2
ref_strains=$3
currentPath=`pwd`

samples=`ls $data_path`

if [ ! -d "variantsAndProp" ]
then
	mkdir "variantsAndProp"
fi

for samp in $samples
do
	bash sample_getVarAndProp.sh $samp $data_path$samp"/" $currentPath"/variantsAndProp/" $lociDb_path $currentPath"/getVariantAndProp.py"
done

if [ ! -d "strainsAndProp" ]
then
	mkdir "strainsAndProp"
fi

echo ""
echo "*********** Solving ILP to find strains and their proportions in each sample ***************"
echo ""
python getStrAndProp.py --data $currentPath"/variantsAndProp" --ref $ref_strains --output $currentPath"/strainsAndProp"

echo ""
echo "Script done."
echo "Created folder variantsAndProp which contains variants identified and their proportions for each sample"
echo "Created folder strainsAndProp which contains strains and their proportions for each sample"


