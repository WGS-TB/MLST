# Illuminating the diversity of Borrelia Burgdorferi in tick samples
All required scripts and files are in pipeline folder (except for sample reads as the files are big). Download the reads for each sample using command `fastq-dump --split-3 [--sample number]`, then place it in a folder named by the sample name and place this folder under `data`.
Once data are ready, run the following commands in the pipeline folder:
```
bash borreliaPipeline.sh [--path to data folder] [--path to loci database folder] [--path to reference strains file]
```
For example,
```
bash borreliaPipeline.sh ~/Documents/Borrelia/pipeline/data/ ~/Documents/Borrelia/pipeline/loci_db/ ~/Documents/Borrelia/pipeline/strain_ref.txt
```
