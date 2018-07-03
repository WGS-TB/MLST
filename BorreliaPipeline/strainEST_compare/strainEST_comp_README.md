# How to run simulations described in our paper?
We assume that all bioinformatics tools that we described in our main page of the github repo is in the user's bin file. We ran all of the simulations on a cluster, with the maximum resource used being 32 cores, 32GB RAM, walltime 6 hours. After and while running the simulations, lots of intermediate files are produced in the current directory for reproducibility and sanity checking. 

1. 
```
Run_StrainEST.py [-h] [-n NUMOFITER] [-st MODTYPE] [-Ed EDITDIST]
                        [-nm NUMMUT] [-em2v EM2V] [-ns NS]
```

`-n` is the number of iterations, `-st` is the EvoMod type, `-Ed` is the maximum edit distance for mutation, `-nm` is the number of mutation event for allele mutation, `em2v` describes EM2e or EM2n. `-ns` is the parameter for benchmark simulation where only existing strains are simulated, it refers to the number of strains.

2. To run benchmark simulation:
* For k=1 (number of existing strains to simulate), `python Run_StrainEST.py -nm 2 -ns 1`.
* For k=2, `python Run_StrainEST.py -nm 2 -ns 2`.

3. To run evoMods simulation:
* EvoMod1: `python Run_StrainEST.py -st EvoMod_1 -nm 2`
* EvoMod2: `python Run_StrainEST.py -st EvoMod_2 -nm 2`
* EvoMod2e and EvoMod2n: `python Run_StrainEST.py -st EvoMOd_2 -nm 2 -em2v [existing/new]` (input accordingly for `em2v`. ** Remember to change the seed to 1992 in Run_StrainEST.py when running EvoMod2e and EvoMod2n **

All simulations will produce `SDP_Output` and `Results` folders, where `SDP_Output` contains some intermediate files corresponding to the output of our algorithm and `Results` contain the statistics.
