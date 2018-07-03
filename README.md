# Deconvoluting the Diversity of Within-host Pathogen Strain in a MLST Framework
## Pre-requisites
1) We are using Python 2.7. You will only need the folder **BorreliaPipeline** to run the pipeline. This can be done in the following `git` commands:
```
git init
git remote add origin https://github.com/WGS-TB/MLST
git fetch origin
git checkout origin/master -- BorreliaPipeline
```

2) SRA toolkit: To download sample reads. Please download the toolkit here: https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/. You will need `fastq-dump` to download the sample reads. Our samples are downloaded using `fastq-dump` v2.8.2. Copy `fastq-dump` into your bin folder (ours: /usr/local/bin) or provide the path to `fastq-dump` to `download_samples.py` script. The python script calls the `fastq-dump` command to download sample reads. 

3) Bowtie v0.12.7: For read mapping. Please download it here: https://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/. You will need both `bowtie` and `bowtie-build`, copy these into your bin folder or provide the path to the folder containing `bowtie` and `bowtie-build` to `mapSamples.py` script. The python script calls `bowtie` and `bowtie-build` commands to map reads.

4) CPLEX Python API: For solving ILP. You will need academic license to download CPLEX. Once you installed CPLEX, follow the instructions here to install CPLEX-Python modules: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html 

5) ART: For simulating reads. Please download it here: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/. You will need `art_illumina`, copy this into your bin file or provide the path for `art_illumina` to `run_sim.py`. This script runs the allele simulation.

6) Samtools v1.3.1: For extracting information from sam files. Please download it here: https://sourceforge.net/projects/samtools/files/samtools/1.3.1/. You will need `samtools`, copy it into your bin folder or provide the path to `samtools` to `mapSamples.py`. The python script calls the `samtools` command while extracting information from sam files after mapping. 

7) In computing edit distance, we used the package by https://github.com/aflc/editdistance and we hereby acknowledge Hiroyuki Tanaka, where the algorithm is proposed by Heikki Hyyrö, "Explaining and extending the bit-parallel approximate string matching algorithm of Myers", (2001). Possible Python packages you will need: numpy, pandas, matplotlib, linecache, argparse . Please use `pip install` to install these packages if you do not have any of the packages.

8) We are not sure if we cover all the Python packages here, please let us know if we miss any!

## Instructions to run the pipeline
1) The required scripts and files are in `pipeline` folder in `BorreliaPipeline` folder (except for sample reads as the files are big). 

2) Run `python download_samples.py [-h] [-c CMD]` in the `pipeline` folder to download the samples' reads (downloaded samples are based on the SRR_Acc_List.txt). `download_samples.py` has an optional argument `-c`, which sets the path for `fastq-dump` command needed to download the reads for Borrelia samples. By default, it assumes `fastq-dump` is in your bin folder. After running this script, it will create a `data` folder which contains folders for each sample.

3) Since we are using Bowtie to map the sample reads, we run `python mapSamples.py [-h] [-c NUMOFCORES] [-b BOWTIE] [-s SAMTOOLS]` to map the sample reads. `-c` specifies the number of core for running Bowtie, default is 4. `-b` specifies the path to the folder containing two commands `bowtie` and `bowtie-build`, default assumes these two commands are in the user's bin folder. `-s` specifies the path of `samtools`, default assumes the command is in user's bin folder. It will be easier to understand with a working example: `python mapSamples.py -c 2 -b ../bowtieFolder/ -s ../../samtools`, for example. The script will create a `variantsAndProp` folder containing a folder for each sample, where each of these folder contains a list of `sampleX_geneY_paired_reads.txt` text files, in which the text file contains information needed for our first stage ILP pipeline. These information are extracted from `.sam` file.

4) To solve the allele diversity problem, run the following script in the pipeline folder:
```
python alleleDiversity.py [-h] [-s SAMPLE]
```
If you want to run on a specific sample, you can indicate using `-s`. Otherwise it will run on all samples in `variantsAndProp`. For example, you could run `python alleleDiversity.py -s SRR2034333`. This script will output some intermediate results to the screen and create 8 CSV files (for each locus) in each of the sample folder in `variantsAndProp`, named as `geneX_proportions.csv` which contains the alleles and their proportions identified at `geneX`.

The script will output some intermediate results for you to keep track.

5) To solve the global strain diversity problem, run the following script:
```
python globalStrainDiversity.py [-h] [-o OUTPUT] [-oc OBJECTIVECOMPONENT] [-timelim TIMELIMIT] [-g GAP] [-pathToDistMat PTDM]
```
Specify `-o`, the folder name for the results to be stored. Default will be `strainsAndProp`. The `-oc` argument specifies different formulation of the objective function, 'all' includes all components of the objective function, 'noPropAndErr' omits the proportion and error component. `-g` specifies the relative gap in % to stop the solver, default is 5%. `-timelim` specifies the maximum time limit in seconds before stopping the solver, default is 600 seconds. The program will stop if and only if `-g` and `-timelim` are both satisfied. For example, if `-g` is 5 and `timelim` is 600, if the relative gap is above 5% and the 600 seconds have passed, the program will not stop until gap of 5% is reached. `-pathToDistMat` is the folder name containing the distance matrices for each gene, which is provided and named 'editDist'.

The script will output some intermediate results for you to keep track. It will create a new folder in your current directory with name specified by `-o`. The folder contains a list of `sampleX_strainsAndProportions.csv` with sampleX corresponds a particular sample, and the csv file contains the strains and their proportions identified in sampleX.

## Instructions to run the allele diversity ADP simulation
1) The required scripts and files to run the simulation are in the `alleleSimulation` folder in the `BorreliaPipeline` folder.

2) In the `alleleSimulation` folder, you only have to run 
```
python run_sim.py [-h] [-i NUMOFITER] [-f SIMULATIONRESULTFOLDER] [-c COVERAGE] [-e EDITDIST] [-fp FULLPAPER] [-b BOWTIE] [-s SAMTOOLS] [-a ART]
``` 
`-i` specifies the number of simulations on each gene, default is 40. `-f` specifies the name of the folder to store the results, default is simulation_results. `-c` specifies the coverage to test on, default is 30. `-b` specifies the path to the folder containing `bowtie` and `bowtie-build` commands, default assumes that both commands are in user's bin folder. `-s` specifies the path to `samtools`, default assumes it's in your bin folder. `-a` specifies the path to `art_illumina`, default assumes that it's in the user's bin folder. `-fp` set to True means we run a 2D parameters experiments(coverage and editDist), which is presented in our paper. In the results folder, it will contain .csv and .png files representing different statistics. 

3) You will need to install `kallisto` and the link is provided here: https://pachterlab.github.io/kallisto/download. To run simulation using kallisto, run
```
run_kallisto_sim.py [-h] [-i NUMOFITER] [-f SIMULATIONRESULTFOLDER] [-c COVERAGE] [-a ART] [-k KAL] [-fp FULLPAPER]                                                     
```

Instructions are similar to `run_sim.py`, with`-k` refers to the path to `kallisto` and default assumes it is in user's bin folder.

4) To run novel allele experiment, run
```
python run_loo_sim.py [-h] [-i NUMOFITER] [-f SIMULATIONRESULTFOLDER] [-c COVERAGE] [-b BOWTIE] [-s SAMTOOLS] [-a ART]
```

Instructions are the same as running `run_sim.py`

## Instructions to run the SDP and full pipeline simulation
1) In the `BorreliaPipeline` folder, `strainAndFullPipeline_sim` contains the scripts and files to run the SDP simulation(with no reads involved) and full pipeline simulation for our method. `strainEST_compare` contains the scripts to run the full pipeline simulation which compares with strainEST. The link to install strainEST is here: https://github.com/compmetagen/strainest. Details about running these simulations are written in separate README files instead in each of the folder.

## For Minimum Spanning Tree Analysis
A sample script to compute the MST is in the `pipeline` folder and can be altered for your own analysis.
