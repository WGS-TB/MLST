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
#samples = ["SRR203433{}".format(i) for i in range(3,7)]
#samples = ["SRR2034333", "SRR2034334", "SRR2034335"]
ap = argparse.ArgumentParser()
ap.add_argument("-s", "--sample", required = True)
args = vars(ap.parse_args())
samples = [args["sample"]]

#currentpath= /pipeline/variantsAndProp
os.chdir("variantsAndProp")

''' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Predicting Variant and Proportions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$''' 
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
        print("..... Predicting variants and computing their proportions .....")
        print("")
        solutionsAndProp_dict = gvp.getVarAndProp(gene,"{0}_{1}_paired_reads.txt".format(samp, gene), samp )
        gene_solProp_dict[gene] = solutionsAndProp_dict
    
        print("Variants and proportions: \n{}".format(solutionsAndProp_dict[0]))
        #write proportions to file
        with open(gene+'_proportions.csv', "wb") as writeFile:
            writer = csv.writer(writeFile)
            for (key, val) in (solutionsAndProp_dict[0]).items():
                writer.writerow([key, val])
                
    #currentpath=/pipeline/variantsAndProp
    os.chdir("..")

print("")
print("Script done.")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
    
