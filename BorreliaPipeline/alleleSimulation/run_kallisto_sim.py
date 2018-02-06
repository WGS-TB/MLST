#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 13:15:44 2017

@author: glgan


"""

import os
import sys
import argparse
#import simulation_functions as sim
import Kallisto_Sim as sim
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
matplotlib.style.use('ggplot')


def construct_boxPlot(csv,type_of_data,name,coverage, iteration):
    if type_of_data == 'Recall':
        label = type_of_data
        limit = [-0.1,1.1]
    elif  type_of_data == 'Precision':
        label = type_of_data
        limit = [-0.1,1.1]
    elif type_of_data == 'DiffObjVal':
        label="Difference in Objective Value: Predicted - True"
        limit = [-10, 10]
    else:
        label = 'Total Variation Distance'
        limit = [-1,101]
    #read in the dataframe
    df = pd.read_csv(csv, sep='\t')
    df.drop(df.columns[[0]], axis=1,inplace=True)
    #create a plot object
    plt.figure()
    #plot the data
    ax = df.plot(kind='box',xlim=[0,9],ylim=limit,title="{0}X Coverage Simulation: {1} simulations each gene".format(coverage, iteration), grid=True)
    ax.set_xlabel('Loci')
    ax.set_ylabel(label)
    plt.xticks([1,2,3,4,5,6,7,8],['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA'])
    plt.savefig(name)
    plt.close('all')

def writeToCsv(fileName, simulationFolderPath, outputFolderName, dataList):
    dataDF = pd.DataFrame(dataList)
    dataDF = dataDF.T
    
    dataDF.to_csv("{0}/{1}/{2}".format(simulationFolderPath, outputFolderName, fileName), sep='\t')
    
start = time.time()
genes = sorted(os.listdir("sim_data"))
genes = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
originalPath=os.getcwd()
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--numOfIter", required = False, default = 40, type=int, help="Number of simulations for each gene")
ap.add_argument("-f", "--simulationResultFolder", required=False, default="simulation_results", help="Folder name to store simulation results, folder with this name will be created in current directory. Default creates a folder named `simulation_results` in current directory.")
ap.add_argument("-c", "--coverage", required=False, default=30,type=int, help="Coverage to test on. Default is 30")
#ap.add_argument("-b", "--bowtie", required=False, default="",help="Path to folder containing bowtie and bowtie-build. Default assumes both bowtie and bowtie-build in user's bin file")
#ap.add_argument("-s", "--samtools", required=False, default="samtools",help="Path to samtools. Default assumes in user's bin file.")
ap.add_argument("-a", "--art", required=False, default="art_illumina",help="Path to art_illumina. Default assumes in user's bin file")

#ap.add_argument("-p", "--proportionMethod", required=True)
args = vars(ap.parse_args())

sim_folder = os.path.abspath(args["simulationResultFolder"])
if args["bowtie"] == "":
    bt = ""
else:
    bt=os.path.abspath(args["bowtie"]) +"/"
    
if args["samtools"] != "samtools":
    samTools = os.path.abspath(args["samtools"])
else:
    samTools = args["samtools"]
    
if args["art"] != "art_illumina":
    art = os.path.abspath(args["art"])
else:
    art = args["art"]

#Make directory for simulation results
directoriesHere = [d for d in os.listdir(".") if os.path.isdir(d)]
if not os.path.isdir(sim_folder):
    os.mkdir(sim_folder)
    
#Run the simulation
os.chdir("sim_data")
originalSTDOut = sys.stdout

#Output statistics as csv file
precision_list = list()
recall_list = list()
totalVarDist_count_list = list()
totalVarDist_bayes_list = list()
for locus in genes:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulating locus {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(locus))
    os.chdir(locus)
    currentDir = os.listdir(".")
    readsFiles = [i for i in currentDir if i.endswith("reads.txt")]
    faFiles = [i for i in currentDir if i.endswith(".fa")]
    
    #Remove previous simulation files if exist
    if len(readsFiles) != 0:
        for i in readsFiles:
            os.remove(i)
        print("==== Removed previous reads.txt files ====")
            
    if len(faFiles) != 0:
        for i in faFiles:
            os.remove(i)
        print("==== Removed previous .fa files ====")
    
    #Function to run simulation imported    
    #precision, recall, diff_obj_vals, totalVarDist_count= sim.simulation(locus,args["numOfIter"],originalPath, args["simulationResultFolder"], args["coverage"],bt,samTools,art)
    precision, recall, totalVarDist_count= sim.simulation(locus,args["numOfIter"],originalPath, args["simulationResultFolder"], args["coverage"])
    precision_list.append(precision)
    recall_list.append(recall)
    totalVarDist_count_list.append(totalVarDist_count)
#    totalVarDist_bayes_list.append(totalVarDist_bayes)
    
    sys.stdout = originalSTDOut

    os.chdir("..")
    print("")
    
#Write these statiscs as csv file and plot it
writeToCsv("{0}X_precision.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], precision_list)
construct_boxPlot("{0}/{1}/{2}X_precision.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "Precision", "{0}/{1}/{2}X_precision.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )
writeToCsv("{0}X_recall.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], recall_list)
construct_boxPlot("{0}/{1}/{2}X_recall.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "Recall", "{0}/{1}/{2}X_recall.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )
writeToCsv("{0}X_totalVarDist_count.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], totalVarDist_count_list)
construct_boxPlot("{0}/{1}/{2}X_totalVarDist_count.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "TotalVarDist", "{0}/{1}/{2}X_totalVarDist_count.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )
#writeToCsv("{0}X_totalVarDist_bayes.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], totalVarDist_bayes_list)
#construct_boxPlot("{0}/{1}/{2}X_totalVarDist_bayes.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "TotalVarDist", "{0}/{1}/{2}X_totalVarDist_bayes.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )

#Summarize the results for all genes
os.chdir("{0}/{1}/".format(originalPath, args["simulationResultFolder"]))
if "allGenes_summary_stats.txt" in os.listdir("."):
    os.remove("allGenes_summary_stats.txt")

geneOutputTxt = [f for f in os.listdir(".") if f.endswith("output_stats.txt")]

for file in geneOutputTxt:
    summarize_cmd = "sed '/SUMMARY STATISTICS*/,$!d' {0} >> allGenes_summary_stats.txt".format(file)    
    os.system(summarize_cmd)

os.chdir(originalPath)  
print("Time taken: {} min".format( (time.time() - start)/60 ))
