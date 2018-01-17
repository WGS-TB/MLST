#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 13:02:22 2017

@author: glgan
"""

import os
import pandas as pd
import time
import argparse

start_time = time.time()

#currentpath = /pipeline/
currentPath = os.getcwd()
lociDb_path = currentPath + "/loci_db/"
ref_strains = currentPath + "/strain_ref.txt"
loci = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]
reference = pd.read_csv(currentPath+"/strain_ref.txt",sep="\t",usecols=range(1,len(loci)+1))
for name in loci:
    reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)

ap = argparse.ArgumentParser()
ap.add_argument("-d", "--data", required=True)
ap.add_argument("-s", "--sample", default='predAll', required = False)
ap.add_argument("-c", "--numOfCores", default=4, required=False)
args = vars(ap.parse_args())

#If particular sample is specified, run on that only. Otherwise, run on all
if args["sample"] != "predAll":
    samples = [args["sample"]]
else:
    samples = os.listdir(os.path.abspath(args["data"]))

#currentpath= /pipeline/variantsAndProp
if not os.path.exists("variantsAndProp"):
    os.mkdir("variantsAndProp")
os.chdir("variantsAndProp")

''' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Predicting Variant and Proportions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$''' 
for samp in samples:
    #/pipeline/variantsAndProp/samp
    if not os.path.exists(samp):
        os.mkdir(samp)
    
    os.chdir(samp)
    print("")
    print("===========================================================================================")
    print("============================= Now processing {} ==================================".format(samp))
    print("===========================================================================================")
    
    for gene in loci:
        print("")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("")
        
        print("..... Predicting variants and computing their proportions .....")
        print("")
        
        firstRead = "{0}/{1}/{1}_1.fastq".format(args["data"], samp)
        secondRead = "{0}/{1}/{1}_2.fastq".format(args["data"], samp)
        kal_cmd = "kallisto quant -t {0} -i {1}{2}.idx -o {2} {3} {4}".format(args["numOfCores"], lociDb_path, gene, firstRead, secondRead)
        os.system(kal_cmd)
        
        #Summarize the proportions based on abundace.tsv file
        abundanceDF = pd.read_csv("{}/abundance.tsv".format(gene), sep='\t')
        varAndPropDF = abundanceDF.loc[:, ['target_id', 'est_counts']]
        varAndPropDF.rename(columns={'target_id':'Allele', 'est_counts':'Proportion'}, inplace=True)
        varAndPropDF = varAndPropDF[varAndPropDF['Proportion'] != 0.0]
        varAndPropDF['Proportion'] = (varAndPropDF["Proportion"]/varAndPropDF["Proportion"].sum())
        varAndPropDF = varAndPropDF.reset_index(drop=True)
        varAndPropDF = varAndPropDF[varAndPropDF.loc[:, "Proportion"] > 0.01]
        varAndPropDF['Proportion'] = (varAndPropDF["Proportion"]/varAndPropDF["Proportion"].sum())
        varAndPropDF.to_csv("{0}_proportions.csv".format(gene), header=False, index=False)
        
    #/pipeline/variantsAndProp
    os.chdir("..")
        
#currentpath=/pipeline/
os.chdir("..")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
    
