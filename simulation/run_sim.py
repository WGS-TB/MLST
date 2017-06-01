#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 13:15:44 2017

@author: glgan
"""

import os
import sys
import argparse
import simulation_functions as sim

genes = os.listdir("sim_data")
originalPath=os.getcwd()
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--numOfIter", required = False, default = 40, type=int)
ap.add_argument("-f", "--simulationResultFolder", required=False, default="simulation_results")
ap.add_argument("-c", "--coverage", required=False, default=30,type=int)
ap.add_argument("-p", "--proportionMethod", required=True)
args = vars(ap.parse_args())

#Make directory for simulation results
directoriesHere = [d for d in os.listdir(".") if os.path.isdir(d)]
if args["simulationResultFolder"] not in directoriesHere:
    os.mkdir(args["simulationResultFolder"])

#Run the simulation
os.chdir("sim_data")
originalSTDOut = sys.stdout
for locus in genes:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulating locus {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(locus))
    os.chdir(locus)
    currentDir = os.listdir(".")
    outFiles = [i for i in currentDir if i.endswith(".out")]
    readsFiles = [i for i in currentDir if i.endswith("reads.txt")]
    faFiles = [i for i in currentDir if i.endswith(".fa")]
    
    #Remove previous simulation files if exist
    if len(outFiles) != 0:
        for i in outFiles:
            os.remove(i)
        print("==== Removed previous .out files ====")
            
    if len(readsFiles) != 0:
        for i in readsFiles:
            os.remove(i)
        print("==== Removed previous reads.txt files ====")
            
    if len(faFiles) != 0:
        for i in faFiles:
            os.remove(i)
        print("==== Removed previous .fa files ====")
    
    #Function to run simulation imported    
    sim.simulation(locus,args["numOfIter"],originalPath, args["simulationResultFolder"], args["coverage"], args["proportionMethod"])
    sys.stdout = originalSTDOut
    
    os.chdir("..")
    print("")

#Summarize the results for all genes
os.chdir("{0}/{1}/".format(originalPath, args["simulationResultFolder"]))
if "allGenes_summary_stats.txt" in os.listdir("."):
    os.remove("allGenes_summary_stats.txt")

geneOutputTxt = [f for f in os.listdir(".") if f.endswith("output_stats.txt")]

for file in geneOutputTxt:
    summarize_cmd = "sed '/SUMMARY STATISTICS*/,$!d' {0} >> allGenes_summary_stats.txt".format(file)    
    os.system(summarize_cmd)

os.chdir(originalPath)    