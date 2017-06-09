#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 13:02:22 2017

@author: glgan
"""

import os
import getVariantAndProp as gvp
import getStrAndProp as gsp
import numpy as np
import itertools
import pandas as pd
import csv

#currentpath = /pipeline/
currentPath = os.getcwd()
data_path = currentPath +"/data/"
lociDb_path = currentPath + "/loci_db/"
ref_strains = currentPath + "/strain_ref.txt"
loci = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]
reference = pd.read_csv(currentPath+"/strain_ref.txt",sep="\t",usecols=range(1,len(loci)+1))
for name in loci:
    reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)

samples = [i for i in os.listdir(data_path)]
#samples = ["SRR2034333"]

if not os.path.exists("variantsAndProp"):
    os.mkdir("variantsAndProp")

#currentpath= /pipeline/variantsAndProp
os.chdir("variantsAndProp")

''' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4 Predicting Variant and Proportions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$''' 
for samp in samples:
    if not os.path.exists(samp):
        os.mkdir(samp)
    
    #currentpath=/pipeline/variantsAndProp/samp
    os.chdir(samp)
    print("")
    print("===========================================================================================")
    print("============================= Now processing {} =================================".format(samp))
    print("===========================================================================================")
    
    #Key=gene, value=dictionary which contains the information about multiple optimal solution
    gene_solProp_dict = dict()
    for gene in loci:
        print("")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("")
        #building index for gene
        buildIndexCmd = "bowtie-build {0}{1}.fas {1}_bowtie >/dev/null".format(lociDb_path, gene)
        os.system(buildIndexCmd)
        print("..... Done building index .....")
        
        #map sample to locus
        print("..... Mapping sample reads to variants .....")
        print("")
        bowtieMapCmd = "bowtie -a -v 3 -p 8 {0}_bowtie -1 {1}{2}_1.fastq -2 {1}{2}_2.fastq {2}_{0}.out >/dev/null".format(gene, data_path+samp+"/", samp)
        os.system(bowtieMapCmd)
        
        #create reads table
        print("..... Summarizing Bowtie mapping information into reads table .....")
        print("")
        readDotOutFile = open("{0}_{1}.out".format(samp, gene))
        writefile = open("{0}_{1}_reads.txt".format(samp, gene), "w")
        for line in readDotOutFile:
            fields = line.strip("\t").split()
            if len(fields) < 10:
                mm_info = ""
            else:
                mm_info = fields[9]
            writefile.write(fields[0] + "\t" + fields[4] + "\t" + mm_info + "\n")
        readDotOutFile.close()
        writefile.close()
        
        #Generate matrix, predict variants and their proportions
        print("..... Predicting variants and computing their proportions .....")
        print("")
        solutionsAndProp_dict = gvp.getVarAndProp(gene,"{0}_{1}_reads.txt".format(samp, gene) )
        gene_solProp_dict[gene] = solutionsAndProp_dict
    
    #Remove intermediate files
    filesRemove = [f for f in os.listdir(".") if (f.endswith(".ebwt") or f.endswith(".out"))]
    for f in filesRemove:
        os.remove(f)
    
    #Choose an optimal solution which maximizes number of existing strains, similar flavour as 2nd ILP
    print(''' Picking a good set of variants among the multiple optimal solutions ''')
    #gene_keys = [[0,1], [0,1,2],...] indices of solutions in each gene
    gene_keys = [gene_solProp_dict[gene].keys() for gene in loci ]
    #Combination of multiple optimal solutions across all genes
    combinationsTuple = [comb for comb in itertools.product(*gene_keys)]
    objValue_list = list()
    track = 0
    for comb in combinationsTuple:
        print("xxxxxxxxxxxxxxxxx Combination : {} xxxxxxxxxxxxxxxxxxxxxxxxxxxx".format(track))
        comb_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb)}
        objval = gvp.maxExistingStr(samp, loci, comb_dict, reference)
        objValue_list.append(objval)
        track += 1
    
    print("Objective Value: {}".format(objValue_list))
    #Choose the combination which has the lowest objective value
    minObjValIndex_list = np.argwhere(objValue_list == np.amin(objValue_list))
    minObjValIndex_list = minObjValIndex_list.flatten().tolist()
    if len(minObjValIndex_list) > 1:
        print("You have more than 1 solution having same objective value")
    
    minObjValIndex = minObjValIndex_list[0]
    comb_minObjVal = combinationsTuple[minObjValIndex]
    comb_minObjVal_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb_minObjVal)}
    
    #Print and write to file
    for gene in comb_minObjVal_dict.keys():
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("Variants and proportions: \n{}".format(comb_minObjVal_dict[gene]))
        #write proportions to file
        with open(gene+'_proportions.csv', "wb") as writeFile:
            writer = csv.writer(writeFile)
            for (key, val) in (comb_minObjVal_dict[gene]).items():
                writer.writerow([key, val])
    #currentpath=/pipeline/variantsAndProp
    os.chdir("..")


''' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Predict strains and their proportions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'''
#currentpath=/pipeline/
os.chdir("..")
if not os.path.exists("strainsAndProp"):
    os.mkdir("strainsAndProp")

print("")
print("******************** Solving ILP to find strains and their proportions in each sample **************************")
print("")
gsp.strainSolver(currentPath+"/variantsAndProp", currentPath+"/strain_ref.txt", currentPath+"/strainsAndProp", loci)

print("")
print("Script done.")
print("Created folder variantsAndProp which contains variants identified and their proportions for each sample")
print("Created folder strainsAndProp which contains strains and their proportions for each sample")
    
