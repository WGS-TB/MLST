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
import argparse
import cplex

def writePropToCsv(aDict):
    for gene in aDict.keys():
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("Variants and proportions: \n{}".format(aDict[gene]))
        #write proportions to file
        with open(gene+'_proportions.csv', "wb") as writeFile:
            writer = csv.writer(writeFile)
            for (key, val) in (aDict[gene]).items():
                writer.writerow([key, val])
                
#def localILP(samp, aTuple, aDict, loci, reference, option):
#    track = 1
#    objValue_list = list()
#    print("\nNumber of combinations to run: {}\n".format(len(aTuple)))
#    for comb in aTuple:
#        print("\nxxxxxxxxxxxxxxxxx Combination : {} xxxxxxxxxxxxxxxxxxxxxxxxxxxx\n".format(track))
#        comb_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, comb)}
#        objVal = gvp.maxExistingStr(samp, loci, comb_dict, reference, option)
#        objValue_list.append(objVal)
#        track += 1
#        
#    print("Objective Value: {}".format(objValue_list))
#    #Choose the combination which has the lowest objective value
#    minObjValIndex_list = np.argwhere(objValue_list == np.amin(objValue_list))
#    minObjValIndex_list = minObjValIndex_list.flatten().tolist()
#    if len(minObjValIndex_list) > 1:
#        print("@@@@@@@@@@@@@@@@@@@@@@@ You have more than 1 distribution having same objective value @@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
#    
#    minObjValIndex = minObjValIndex_list[0]
#    comb_minObjVal = aTuple[minObjValIndex]
#    comb_minObjVal_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, comb_minObjVal)}
#    
#    return comb_minObjVal_dict

def localMinimizer(samp, aTuple, aDict, loci, reference, option):
    track = 1
    objValue_list = list()
    print("\nNumber of combinations to run: {}\n".format(len(aTuple)))
    for comb in aTuple:
        print("\nxxxxxxxxxxxxxxxxx Combination : {} xxxxxxxxxxxxxxxxxxxxxxxxxxxx\n".format(track))
        comb_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, comb)}
        solution_dict, ilp_objective_dict, data, strains = gvp.localILP(samp, loci, comb_dict, reference)
        feasible_sol = list()
        lp_objective_dict = dict()
#        print solution_dict
        infeasibility = 0
        for i in solution_dict.keys():
            try:
                objvalue = gvp.localLP(solution_dict[i], data, strains, reference, loci)
                feasible_sol.append(i)
                lp_objective_dict[i] = objvalue
            except cplex.exceptions.errors.CplexSolverError as e:
                infeasibility += 1
        
        if infeasibility == len(solution_dict):
            return -1
        
        min_obj = np.inf
        for j in feasible_sol:
            if (ilp_objective_dict[j] + lp_objective_dict[j]) < min_obj:
                min_obj = ilp_objective_dict[j] + lp_objective_dict[j]
        
        objValue_list.append(min_obj)
        print("Objective value: {}".format(min_obj))
        track += 1
        
    print("Objective Value: {}".format(objValue_list))
    #Choose the combination which has the lowest objective value
    minObjValIndex_list = np.argwhere(objValue_list == np.amin(objValue_list))
    minObjValIndex_list = minObjValIndex_list.flatten().tolist()
    if len(minObjValIndex_list) > 1:
        print("@@@@@@@@@@@@@@@@@@@@@@@ You have more than 1 distribution having same objective value @@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    
    minObjValIndex = minObjValIndex_list[0]
    comb_minObjVal = aTuple[minObjValIndex]
    comb_minObjVal_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, comb_minObjVal)}
    
    return comb_minObjVal_dict

def compatibleFilter(aTuple):
    compatible_tuples = list()
    for comb in aTuple:
        comb_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, comb)}
        
        temp_boolean = True
        for allele in loci:
            temp_boolean = temp_boolean & reference[allele].isin(comb_dict[allele].keys())
            if sum(temp_boolean) == 0:
                break
            
        if sum(temp_boolean) != 0:
            compatible_tuples.append(comb)
            
    return compatible_tuples
''' ----------------------------------------------------------------- '''
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

#samples = [i for i in os.listdir(data_path)]
ap = argparse.ArgumentParser()
ap.add_argument("-s", "--sample", required = True)
ap.add_argument("-cl", "--objectiveComponent_local", default="all", required=False)
ap.add_argument("-cg", "--objectiveComponent_global", default="all", required=False)
args = vars(ap.parse_args())
samples = [args["sample"]]

#currentpath= /pipeline/variantsAndProp
os.chdir("variantsAndProp")

''' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Predicting Variant and Proportions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$''' 
numOfMinQualSol = list()
for samp in samples:
    #currentpath=/pipeline/variantsAndProp/samp
    os.chdir(samp)
    print("")
    print("===========================================================================================")
    print("============================= Now processing {} ==================================".format(samp))
    print("===========================================================================================")
    
    #Key=gene, value=dictionary which contains the information about multiple optimal solution
    gene_solProp_dict = dict()
    for gene in loci:
        print("")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("")
        
        #Generate matrix, predict variants and their proportions
        print("..... Predicting variants .....")
        print("")
        solutionsAndProp_dict = gvp.getVarAndProp(gene,"{0}_{1}_paired_reads.txt".format(samp, gene), samp )
        gene_solProp_dict[gene] = solutionsAndProp_dict
        numOfMinQualSol.append(len(solutionsAndProp_dict))
        
    #gene_keys = [[0,1], [0,1,2],...] indices of solutions in each gene
    gene_keys = [gene_solProp_dict[gene].keys() for gene in loci ]
    #Combination of multiple optimal solutions across all genes
    combinationsTuple = [comb for comb in itertools.product(*gene_keys)]
    
    #If each gene has a unique optimal solution
    if all(i == 1 for i in numOfMinQualSol):
    #    Print and write to file
        write_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, combinationsTuple[0])}
        writePropToCsv(write_dict)
    else:   #when there are indistinguishable solutions based on quality score
        #Choose an optimal solution which maximizes number of existing strains, similar flavour as 2nd ILP
        print("\n... There are at least one gene having multiple optimal solutions ...")
        print("\n... Picking a good set of variants among the multiple optimal solutions ...")
    
        compatible_tuples = compatibleFilter(combinationsTuple)
        localMinimizer_dict = dict()
        objValue_list = list()
        track = 1
        
        #Only one compatible distribution, output it
        if len(compatible_tuples) == 1:
            localMinimizer_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, compatible_tuples[0])}
            writePropToCsv(localMinimizer_dict)
        else:
            #Only new strains can describe
            if len(compatible_tuples) == 0:
                localMinimizer_dict = localMinimizer(samp, combinationsTuple, gene_solProp_dict, loci, reference, args["objectiveComponent_local"])
            else:
                localMinimizer_dict = localMinimizer(samp, compatible_tuples, gene_solProp_dict, loci, reference, args["objectiveComponent_local"])
        
            if localMinimizer_dict == -1:
                print("No feasibile solutions")
            else:
                writePropToCsv(localMinimizer_dict)
                
    #currentpath=/pipeline/variantsAndProp
    os.chdir("..")

print("")
print("Script done.")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
    
