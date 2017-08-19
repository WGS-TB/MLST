#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  12 13:02:22 2017

@author: glgan
"""

import os
import time
import argparse

def writeReadTable(path, samp, gene, option):
    readOutFile = open("{0}{1}_{2}_{3}NoHeader.sam".format(path, samp, gene, option))
    writefile = open("{0}_{1}_{2}_reads.txt".format(samp, gene, option), "w")
    for line in readOutFile:
        fields = line.strip("\t").split()
        read = fields[0]
        allele = fields[2]
        quality = fields[10]
#        mm = [i for i in fields if i.startswith("XM")][0]  #bowtie2
        mm = [i for i in fields if i.startswith("NM:i:")][0]   #bowtie
        mm_pos = [j for j in fields if j.startswith("MD")][0]
        
        writefile.write(read + "\t" + allele + "\t" + quality + "\t" + mm + "\t" + mm_pos + '\n')
        
    readOutFile.close()
    writefile.close()

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
    temp_start = time.time()
           
    for gene in loci:
        print("")
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~ Gene {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(gene))
        print("")
        #building index for gene
        buildIndexCmd = "bowtie-build {0}{1}.fas {1}_bowtie >/dev/null".format(lociDb_path, gene)
#        buildIndexCmd = "bowtie2-build {0}{1}.fas {1}_bowtie2 >/dev/null".format(lociDb_path, gene)
        os.system(buildIndexCmd)
        print("..... Done building index .....")
        
        #map sample to locus
        print("..... Mapping sample reads to variants .....")
        print("")
        bowtieMapCmd = "bowtie --best --strata -a -v 3 -p {0} {1}_bowtie -1 {2}{3}_1.fastq -2 {2}{3}_2.fastq -S {2}{3}_{1}.sam >/dev/null".format(args["numOfCores"],gene, data_path+samp+"/", samp)
#        bowtieMapCmd = 'bowtie2 -x {1}_bowtie2 -a -p {0} -1 {2}{3}_1.fastq -2 {2}{3}_2.fastq -S {2}{3}_{1}.sam >/dev/null 2>&1'.format(args["numOfCores"],gene, data_path+samp+"/", samp)        
        os.system(bowtieMapCmd)
        
        #handle sam files for writing reads txt file
        mapped_cmd = "samtools view -h -F4 {0}{1}_{2}.sam > {0}{1}_{2}_mapped.sam".format(data_path+samp+"/", samp, gene)
        paired_cmd = "samtools view -F8 {0}{1}_{2}_mapped.sam > {0}{1}_{2}_pairedNoHeader.sam".format(data_path+samp+"/", samp, gene)
#        singleton_cmd = "samtools view -f8 {0}{1}_{2}_mapped.sam > {0}{1}_{2}_singletonNoHeader.sam".format(data_path+samp+"/", samp, gene)
        
        os.system(mapped_cmd)
        os.system(paired_cmd)
#        os.system(singleton_cmd)
        
        #create reads table
        print("..... Summarizing Bowtie mapping information into reads table .....")
        print("")
        writeReadTable(data_path+samp+"/", samp, gene, "paired")
#        writeReadTable(data_path+samp+"/", samp, gene, "singleton")
        samFileRemove = [f for f in os.listdir(data_path+samp) if f.endswith(".sam")]
        for f in samFileRemove:
            os.remove(data_path+samp+"/"+f)
    
    #Remove intermediate files
    filesRemove = [f for f in os.listdir(".") if (f.endswith(".ebwt") or f.endswith(".bt2") or f.endswith("*.out"))]
    for f in filesRemove:
        os.remove(f)
        
    #currentpath=/pipeline/variantsAndProp
    os.chdir("..")
        
    print("Time take for sample {0}: {1}".format(samp, (time.time() - temp_start)/3600 ))
    
print("Total time taken for Bowtie mapping: {} hr(s)".format((time.time() - start_time)/3600))
    
