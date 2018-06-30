# How to run simulations described in our paper?
We assume that all bioinformatics tools that we described in our main page of the github repo is in the user's bin file

1. 
```
python EvoMod_Experiments.py [-h] [-n NUMOFITER] [-d MASTERDIR] [-Ed EDITDIST] [-nm NUMMUT] [-st MODTYPE] [-et EXPERIMENTTYPE] [-emd EMDDISTMODE] [-art ART]
                             [-pathToDistMat PTDM] [-em2_variant EM2V] 
```

`-n` refers to the number of iterations, `-d` is where the results will be stored, `-Ed` is the maximum edit distance for allele mutation, `-nm` is the number of mutation events, `-st` is the EvoMod type, `-et` is the experiment type, 1=Strain only,2=full pipeline, `-emd` is the distance used for EMD calculations ('snp' or 'hamming'), `-art` is the path to art_illumina tool, `-pathToDistMat` is the name of folder which contains the edit distance matrices(please use `editDist` here to produce our simulation), `em2v` is a variant of EvoMod2 (which are evoMod2e/2n reported in the paper). You may run `python EvoMod_Experiments.py-h` to get the available options and details for each parameters.

2. To run SDP simulation for each evoMod simulation (some default parameters are used and can be changed, such as name of folder):
* EvoMod1: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_1 -et 1 -pathToDistMat editDist`
* EvoMod2: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_2 -et 1 -pathToDistMat editDist`
* EvoMod2e: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_2 -et 1 -pathToDistMat editDist -em2_variant existing`
* EvoMod2e: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_2 -et 1 -pathToDistMat editDist -em2_variant new`
* EvoMod3: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_3 -et 1 -pathToDistMat editDist`

3. To run full pipeline simulation for each evoMod simulation:
* EvoMod1: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_1 -et 2 -pathToDistMat editDist`
* EvoMod2: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_2 -et 2 -pathToDistMat editDist`
* EvoMod2e: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_2 -et 2 -pathToDistMat editDist -em2_variant existing`
* EvoMod2e: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_2 -et 2 -pathToDistMat editDist -em2_variant new`
* EvoMod3: `python EvoMod_Experiments.py -Ed 15 -nm 2 -st Type_3 -et 2 -pathToDistMat editDist`

