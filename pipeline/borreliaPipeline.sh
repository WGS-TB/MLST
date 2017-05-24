#!/bin/bash
#Run this script by 'bash borreliaPipelin.sh' and it will create 2 folders, variantsAndProp and strainsAndProp
currentPath=`pwd`
data_path=$currentPath"/data/"
lociDb_path=$currentPath"/loci_db/"
ref_strains=$currentPath"/strain_ref.txt"

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


