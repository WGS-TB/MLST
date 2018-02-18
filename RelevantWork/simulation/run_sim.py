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

def fullPaper_writeToCsv( stats, simulationFolderPath, outputFolderName, gene_covEd_stats):
    genes = gene_covEd_stats.keys()

    for g in genes:
        df = pd.DataFrame.from_dict(gene_covEd_stats[g], orient='index')
        df.to_csv(os.path.join(simulationFolderPath, outputFolderName, "fullPaper_{0}_{1}.csv".format(g, stats)), sep="\t")    

def fullPaper_avg_writeToCsv(stats, simulationFolderPath, outputFolderName, covEd_stats):
    df = pd.DataFrame.from_dict(covEd_stats, orient='index')
    df.to_csv(os.path.join(simulationFolderPath, outputFolderName, "fullPaper_average_{0}.csv".format(stats)), sep="\t")    

def fullPaper_boxplot(originalPath, sim_result_folder, csv, type_of_data, gene, iteration):
    df = pd.read_csv(os.path.join(originalPath, sim_result_folder, csv), sep="\t")
    df.rename(columns={"Unnamed: 0" : "CovEd"}, inplace=True)

    if type_of_data == 'Recall':
        label = type_of_data
        limit = [-0.1,1.1]
    elif  type_of_data == 'Precision':
        label = type_of_data
        limit = [-0.1,1.1]
    elif type_of_data == 'NumOfOpt':
        label="Number of optimal solutions"
        maxOpt = df.iloc[:, 1:].max().max()
        limit = [0,maxOpt+1]
    else:
        label = 'Total Variation Distance'
        limit = [-1,101]

    #create a plot object
    plt.figure()
    #plot the data
    ax = df.plot(kind='box',xlim=[0,df.shape[0]+1],ylim=limit,title="Locus {}".format(gene), grid=True)
    ax.set_xlabel('Experiments: {} simulations each'.format(iteration))
    ax.set_ylabel(label)
    plt.xticks(range(1, df.shape[0]+1), df["CovEd"].tolist())
    plt.savefig(os.path.join(originalPath, sim_result_folder, "{0}_{1}.png".format(gene, type_of_data)))
    plt.close('all')

def main_simulationFunc(locus, iteration, originalPath, simulationResultFolder, cov, ed, bt, samtools, art):
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
    precision, recall, totalVarDist_count, numOfOptSol = sim.simulation(locus,iteration,originalPath, simulationResultFolder, cov, ed ,bt,samTools,art)
    precision_list.append(precision)
    recall_list.append(recall)
    totalVarDist_count_list.append(totalVarDist_count)
    
    sys.stdout = originalSTDOut

    os.chdir("..")
    print("")

    return precision, recall, totalVarDist_count, numOfOptSol
    
start = time.time()
genes = sorted(os.listdir("sim_data"))
genes = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
originalPath=os.getcwd()
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--numOfIter", required = False, default = 40, type=int, help="Number of simulations for each gene")
ap.add_argument("-f", "--simulationResultFolder", required=False, default="simulation_results", help="Folder name to store simulation results, folder with this name will be created in current directory. Default creates a folder named `simulation_results` in current directory.")
ap.add_argument("-c", "--coverage", required=False, default=30,type=int, help="Coverage to test on. Default is 30")
ap.add_argument("-e", "--editDist", required=False, default=10,type=int, help="Maximum edit distance from seed allele to filter alleles for simulation. Default is 5")
ap.add_argument("-fp", "--fullPaper", required=False, default=False, type=bool, help="Full simulation of ADP on 2 dimensional parameters, coverage and maxEditDist. If this is set, then any settings of -e and -c are ignored")
ap.add_argument("-b", "--bowtie", required=False, default="",help="Path to folder containing bowtie and bowtie-build. Default assumes both bowtie and bowtie-build in user's bin file")
ap.add_argument("-s", "--samtools", required=False, default="samtools",help="Path to samtools. Default assumes in user's bin file.")
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
#diffObjVal_list = list()
totalVarDist_count_list = list()
numOfOptSol = list()

#2D parameters
if args["fullPaper"] == True:
    maxEditDist = [5,10,15,20,25]
    maxEditDist=[5,10]
    coverage=[30,100,300]
    coverage=[30]
    gene_covEd_precision = dict()   #key=gene name, value=dictionary, where key=(coverage, editDistance) and value=a list of values
    gene_covEd_recall = dict()
    gene_covEd_tvd = dict()
    gene_covEd_numOfOpt = dict()

for locus in genes:
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulating locus {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(locus))
    os.chdir(locus)
    
    if args["fullPaper"] == False:
        precision, recall, totalVarDist_count, numOfOpt = main_simulationFunc(locus,args["numOfIter"],originalPath, args["simulationResultFolder"], args["coverage"], args["editDist"],bt,samTools,art)
        precision_list.append(precision)
        recall_list.append(recall)
        totalVarDist_count_list.append(totalVarDist_count)
    else:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulating 2D parameters: coverage and edit distances ~~~~~~~~~~~~~~~~~~~~~~~~")
        covEd_prec = dict()
        covEd_rec = dict()
        covEd_tvd = dict()
        covEd_opt = dict()
        for cov in coverage:
            print cov
            for ed in maxEditDist:
                print ed
                precision, recall, totalVarDist_count, numOfOptSol = main_simulationFunc(locus,args["numOfIter"],originalPath, args["simulationResultFolder"], cov, ed,bt,samTools,art)
                covEd_prec[("{}X".format(cov), "{}editDist".format(ed))] = precision
                covEd_rec[("{}X".format(cov), "{}editDist".format(ed))] = recall
                covEd_tvd[("{}X".format(cov), "{}editDist".format(ed))] = totalVarDist_count
                covEd_opt[("{}X".format(cov), "{}editDist".format(ed))] = numOfOptSol
                os.chdir(locus)

        gene_covEd_precision[locus] = covEd_prec        
        gene_covEd_recall[locus] = covEd_rec        
        gene_covEd_tvd[locus] = covEd_tvd        
        gene_covEd_numOfOpt[locus] = covEd_opt
        os.chdir("..")        
    
if args["fullPaper"] == False:    
    #Write these statiscs as csv file and plot it
    writeToCsv("{0}X_precision.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], precision_list)
    #construct_boxPlot("{0}/{1}/{2}X_precision.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "Precision", "{0}/{1}/{2}X_precision.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )
    writeToCsv("{0}X_recall.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], recall_list)
    #construct_boxPlot("{0}/{1}/{2}X_recall.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "Recall", "{0}/{1}/{2}X_recall.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )
    #writeToCsv("{0}X_diffObjVal.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], diffObjVal_list)
    #construct_boxPlot("{0}/{1}/{2}X_diffObjVal.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "DiffObjVal", "{0}/{1}/{2}X_diffObjVal.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )
    writeToCsv("{0}X_totalVarDist_count.csv".format(args["coverage"]), originalPath, args["simulationResultFolder"], totalVarDist_count_list)
    #construct_boxPlot("{0}/{1}/{2}X_totalVarDist_count.csv".format(originalPath, args["simulationResultFolder"], args["coverage"]), "TotalVarDist", "{0}/{1}/{2}X_totalVarDist_count.png".format(originalPath, args["simulationResultFolder"], args["coverage"]), args["coverage"], args["numOfIter"] )

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
else:
    fullPaper_writeToCsv( "precision", originalPath, args["simulationResultFolder"], gene_covEd_precision)
    fullPaper_writeToCsv( "recall", originalPath, args["simulationResultFolder"], gene_covEd_recall)
    fullPaper_writeToCsv( "tvd", originalPath, args["simulationResultFolder"], gene_covEd_tvd)
    fullPaper_writeToCsv( "numOfOpt", originalPath, args["simulationResultFolder"], gene_covEd_numOfOpt)

    for g in genes:
        fullPaper_boxplot(originalPath, args["simulationResultFolder"], "fullPaper_{}_precision.csv".format(g), "Precision", g, args["numOfIter"])    
        fullPaper_boxplot(originalPath, args["simulationResultFolder"], "fullPaper_{}_recall.csv".format(g), "Recall", g, args["numOfIter"])    
        fullPaper_boxplot(originalPath, args["simulationResultFolder"], "fullPaper_{}_tvd.csv".format(g), "Total Variation Distance", g, args["numOfIter"])    
        fullPaper_boxplot(originalPath, args["simulationResultFolder"], "fullPaper_{}_numOfOpt.csv".format(g), "NumOfOpt", g, args["numOfIter"])    
    
