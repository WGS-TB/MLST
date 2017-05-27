# Illuminating the diversity of Borrelia Burgdorferi in tick samples
## Pre-requisites
1) SRA toolkit: To download samples' reads. Please download the toolkit here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/. You will need `fastq-dump` to download the reads. Copy `fastq-dump` into your bin file.

2) Bowtie: For read mapping. Please download it here: http://bowtie-bio.sourceforge.net/index.shtml. You will need both `bowtie` and `bowtie-build`, copy these into your bin file.

3) CPLEX Python API: For solving ILP. You will need academic license to download CPLEX. Once you installed CPLEX, follow the instructions here to install CPLEX-Python modules: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html 

4) Python packages such as numpy and pandas. Please use `pip install` to install these packages.

##Instructions to run the pipeline
1) All required scripts and files are in pipeline folder (except for sample reads as the files are big). 

2) Run `python download_samples.py` to download the samples' reads (downloaded samples are based on the SRR_Acc_List.txt)

3) Once data are ready, run the following commands in the pipeline folder:
```
bash borreliaPipeline.sh
```

4) After the script is done, two folders will be created i.e. variantsAndProp and strainsAndProp. The variantsAndProp folder contains the variants identified and their proportions for each sample. The strainsAndProp folder contains the strains identified and their proportions for each sample
