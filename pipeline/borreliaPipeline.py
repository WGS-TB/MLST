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
import time

start_time = time.time()

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
#samples = ["SRR203433{}".format(i) for i in range(3,7)]
#samples = ["SRR2034333", "SRR2034334", "SRR2034335"]

#currentpath= /pipeline/variantsAndProp
os.chdir("variantsAndProp")

''' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Predicting Variant and Proportions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$''' 
for samp in samples:
    #currentpath=/pipeline/variantsAndProp/samp
    os.chdir(samp)
    print("")
    print("===========================================================================================")
    print("============================= Now processing {} =================================".format(samp))
    print("===========================================================================================")
    
    #Key=gene, value=dictionary which contains the information about multiple optimal solution
    gene_solProp_dict = dict()
    minVar_list = list()
    for gene in loci:
        print("")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("")
        
        #Generate matrix, predict variants and their proportions
        print("..... Predicting variants and computing their proportions .....")
        print("")
        solutionsAndProp_dict, num_minVar_solutions = gvp.getVarAndProp(gene,"{0}_{1}_reads.txt".format(samp, gene), samp )
        print(solutionsAndProp_dict)
        gene_solProp_dict[gene] = solutionsAndProp_dict
        minVar_list.append(num_minVar_solutions)
    
    #Choose an optimal solution which maximizes number of existing strains, similar flavour as 2nd ILP
    print(''' Picking a good set of variants among the multiple optimal solutions ''')
    #gene_keys = [[0,1], [0,1,2],...] indices of solutions in each gene
    gene_keys = [gene_solProp_dict[gene].keys() for gene in loci ]
    #Combination of multiple optimal solutions across all genes
    combinationsTuple = [comb for comb in itertools.product(*gene_keys)]
#    objValue_list = list()
    track = 1
#    for comb in combinationsTuple:
#        print("xxxxxxxxxxxxxxxxx Combination : {} xxxxxxxxxxxxxxxxxxxxxxxxxxxx".format(track))
#        comb_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb)}
#        objval = gvp.maxExistingStr(samp, loci, comb_dict, reference)
#        objValue_list.append(objval)
#        track += 1
#    score_list = list()
    print("Number of combinations using minimum variants only: {}".format(np.prod(minVar_list)))
    print("Number of combinations using minimum variants+likelihood: {}".format(len(combinationsTuple)))
    compatible_tuples = list()
    for comb in combinationsTuple:
        print("xxxxxxxxxxxxxxxxxxxxxxxxx Combination : {} xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx".format(track))
        track += 1
        comb_list = [gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb)]
        comb_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb)}
        
        temp_boolean = True
        for allele in loci:
            temp_boolean = temp_boolean & reference[allele].isin(comb_dict[allele].keys())
            if sum(temp_boolean) == 0:
                break
            
        if sum(temp_boolean) != 0:
            compatible_tuples.append(comb)
    
#    print compatible_tuples
    print("Number of compatible combinations: {}".format(len(compatible_tuples)))
            
#        sol_combinations = list(itertools.product(*comb_list))
#        existing = 0
#        for strain in sol_combinations:
#            existing += sum(reference.isin(strain).all(1))
#        score_list.append(100.0 * (len(sol_combinations) - existing)/(len(sol_combinations)))
#        track += 1
        
#    print("Objective Value: {}".format(objValue_list))
    #Choose the combination which has the lowest objective value
#    minObjValIndex_list = np.argwhere(objValue_list == np.amin(objValue_list))
#    minObjValIndex_list = minObjValIndex_list.flatten().tolist()
#    if len(minObjValIndex_list) > 1:
#        print("You have more than 1 solution having same objective value")
#    
#    minObjValIndex = minObjValIndex_list[0]
#    comb_minObjVal = combinationsTuple[minObjValIndex]
#    comb_minObjVal_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb_minObjVal)}
    
#    print("Score list:{}\n".format(score_list))
#    minScore_list = np.argwhere(score_list == np.amin(score_list))
#    minScore_list = minScore_list.flatten().tolist()
#    if len(minScore_list) > 1:
#        print("You have more than 1 solution having same score")
    
#    minObjValIndex = minScore_list[0]
#    comb_minObjVal = combinationsTuple[minObjValIndex]
#    comb_minObjVal_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb_minObjVal)}
    
    #Print and write to file
#    for gene in comb_minObjVal_dict.keys():
#        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
#        print("Variants and proportions: \n{}".format(comb_minObjVal_dict[gene]))
#        #write proportions to file
#        with open(gene+'_proportions.csv', "wb") as writeFile:
#            writer = csv.writer(writeFile)
#            for (key, val) in (comb_minObjVal_dict[gene]).items():
#                writer.writerow([key, val])
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
#gsp.strainSolver(currentPath+"/variantsAndProp", currentPath+"/strain_ref.txt", currentPath+"/strainsAndProp", loci)

print("")
print("Script done.")
print("Created folder variantsAndProp which contains variants identified and their proportions for each sample")
print("Created folder strainsAndProp which contains strains and their proportions for each sample")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
    
