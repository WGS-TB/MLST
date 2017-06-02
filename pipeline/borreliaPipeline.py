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

#currentpath = /pipeline/
currentPath = os.getcwd()
data_path = currentPath +"/data/"
lociDb_path = currentPath + "/loci_db/"
ref_strains = currentPath + "/strain_ref.txt"
loci = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]

samples = [i for i in os.listdir(data_path)]

if not os.path.exists("variantsAndProp"):
    os.mkdir("variantsAndProp")

#currentpath= /pipeline/variantsAndProp
os.chdir("variantsAndProp")

'''Predicting Variant and Proportions'''
numOfOptSol_list = list()    
for samp in samples:
    if not os.path.exists(samp):
        os.mkdir(samp)
    
    #currentpath=/pipeline/variantsAndProp/samp
    os.chdir(samp)
    print("")
    print("===========================================================================================")
    print("============================= Now processing {} =================================".format(samp))
    print("===========================================================================================")
    
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
        numOfOptSol = gvp.getVarAndProp(gene,"{0}_{1}_reads.txt".format(samp, gene) )
        numOfOptSol_list.append(numOfOptSol)
    
    filesRemove = [f for f in os.listdir(".") if (f.endswith(".ebwt") or f.endswith(".out"))]
    
    for f in filesRemove:
        os.remove(f)
    
    #currentpath=/pipeline/variantsAndProp
    os.chdir("..")

print("Number of ILP needed to solve for 2nd stage: {}".format(np.prod(np.array(numOfOptSol_list))))

'''Predict strains and their proportions'''
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
    
