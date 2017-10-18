#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 13:02:22 2017

@author: glgan
"""

import os
import pipeline_functions as pf
import itertools
import pandas as pd
import csv
import time
import argparse

def writePropToCsv(aDict):
    for gene in aDict.keys():
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("Variants and proportions: \n{}".format(aDict[gene]))
        #write proportions to file
        with open(gene+'_proportions.csv', "wb") as writeFile:
            writer = csv.writer(writeFile)
            for (key, val) in (aDict[gene]).items():
                writer.writerow([key, val])

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
ap.add_argument("-s", "--sample", required = True, help="Sample name")
#ap.add_argument("-lo", "--localOption", required=True, help="Version of optimization program to use. 'mixed': mixed ILP, 'separated': pure ILP + LP")
#ap.add_argument("-timelim", "--timeLimit", required=False, help="Time limit in integer(sec) for cplex solver for mixed ILP. Default: 600sec", default=600)
#ap.add_argument("-g", "--gap", required=False, help="Relative gap tolerance(percent) for cplex solver for mixed ILP. Default: 5", default=5)

#only for MILP
#ap.add_argument("-oc", "--objectiveComponent", required=False, default="all", help="Objective components, only applicable to mixed ILP. Default: 'all'. 'noPropAndErr': Does not include proportion and error in objective function")
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
        solutionsAndProp_dict = pf.getVarAndProp(gene,"{0}_{1}_paired_reads.txt".format(samp, gene), samp )
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
    
        compatible_tuples = pf.compatibleFilter(combinationsTuple, gene_solProp_dict, loci, reference)
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
                localMinimizer_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, combinationsTuple[0])}
            else:
                localMinimizer_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, compatible_tuples[0])}
                
            writePropToCsv(localMinimizer_dict)
#            if len(compatible_tuples) == 0:
#                if args["localOption"] == "mixed":
#                    localMinimizer_dict = pf.localMinimizer(samp, combinationsTuple, gene_solProp_dict, loci, reference, args["objectiveComponent"],args["timeLimit"], args["gap"])
#                else:
#                    localMinimizer_dict = pf.localMinimizer_sep(samp, combinationsTuple, gene_solProp_dict, loci, reference)
#            else: #more than one compatible combinations
#                if args["localOption"] == "mixed":
#                    localMinimizer_dict = pf.localMinimizer(samp, compatible_tuples, gene_solProp_dict, loci, reference, args["objectiveComponent"],args["timeLimit"], args["gap"])
#                else:
#                    localMinimizer_dict = pf.localMinimizer_sep(samp, compatible_tuples, gene_solProp_dict, loci, reference)
            
            #only matters if separated version of local optimization is chosen
#            if localMinimizer_dict == -1:
#                print("No feasibile solutions")
#            else:
#                writePropToCsv(localMinimizer_dict)
                
    #currentpath=/pipeline/variantsAndProp
    os.chdir("..")

print("")
print("Script done.")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
    
