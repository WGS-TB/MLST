# Illuminating the diversity of Borrelia Burgdorferi in tick samples
1) All required scripts and files are in pipeline folder (except for sample reads as the files are big). 

2) Download the reads for each sample using command `fastq-dump --split-3 [--sample number]`, then place it in a folder named by the sample name and place this folder under `data`.

3) Once data are ready, run the following commands in the pipeline folder:
```
bash borreliaPipeline.sh [--path to data folder] [--path to loci database folder] [--path to reference strains file]
```
For example,
```
bash borreliaPipeline.sh ~/Documents/Borrelia/pipeline/data/ ~/Documents/Borrelia/pipeline/loci_db/ ~/Documents/Borrelia/pipeline/strain_ref.txt
```

4) After the script is done, two folders will be created i.e. variantsAndProp and strainsAndProp. The variantsAndProp folder contains the variants identified and their proportions for each sample. The strainsAndProp folder contains the strains identified and their proportions for each sample
