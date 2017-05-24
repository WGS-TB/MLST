# Illuminating the diversity of Borrelia Burgdorferi in tick samples
1) All required scripts and files are in pipeline folder (except for sample reads as the files are big). 

2) You will need SRA toolkit to download the reads. Please download the toolkit here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/. You will need `fastq-dump` to download the reads. Copy `fastq-dump` into your bin file.

3) You will also need Bowtie. Please download it here: http://bowtie-bio.sourceforge.net/index.shtml. You will need both `bowtie` and `bowtie-build`, copy these into your bin file.

4) Download the reads for each sample using command `fastq-dump --split-3 [--sample number]`, then place it in a folder named by the sample name and place this folder under `data`. For example, `/data/SRR2034333/` 

5) Once data are ready, run the following commands in the pipeline folder:
```
bash borreliaPipeline.sh
```

4) After the script is done, two folders will be created i.e. variantsAndProp and strainsAndProp. The variantsAndProp folder contains the variants identified and their proportions for each sample. The strainsAndProp folder contains the strains identified and their proportions for each sample
