#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  12 13:02:22 2017

@author: glgan
"""

import os
import time
import argparse

start_time = time.time()
ap = argparse.ArgumentParser()
ap.add_argument("-c", "--numOfCores", required = False, default = 4, type=int)
args = vars(ap.parse_args())

#currentpath = /pipeline/
currentPath = os.getcwd()
data_path = currentPath +"/data/"
lociDb_path = currentPath + "/loci_db/"
loci = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]

samples = [i for i in os.listdir(data_path)]

#currentpath= /pipeline/variantsAndProp
if not os.path.exists("variantsAndProp"):
    os.mkdir("variantsAndProp")
os.chdir("variantsAndProp")

''' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Mapping Samples $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$''' 
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
        bowtieMapCmd = "bowtie -a -v 3 -p {0} {1}_bowtie -1 {2}{3}_1.fastq -2 {2}{3}_2.fastq {3}_{1}.out >/dev/null".format(args["numOfCores"],gene, data_path+samp+"/", samp)
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
    
    #Remove intermediate files
    filesRemove = [f for f in os.listdir(".") if (f.endswith(".ebwt") or f.endswith(".out"))]
    for f in filesRemove:
        os.remove(f)
        
print("Time taken for Bowtie mapping: {} hr(s)".format((time.time() - start_time)/3600))
    
