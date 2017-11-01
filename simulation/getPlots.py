#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:28:49 2017

@author: glgan
"""

import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import os
matplotlib.style.use('ggplot')

def construct_boxPlot(csv,type_of_data,name,coverage):
    if type_of_data == 'Recall':
        label = type_of_data
        limit = [0,2]
    elif  type_of_data == 'Precision':
        label = type_of_data
        limit = [0,2]
    elif type_of_data == 'DiffObjVal':
        label="Difference in Objective Value: Predicted - True"
        limit = [-10, 10]
    else:
        label = 'Total Variation Distance'
        limit = [-20,100]
    #read in the dataframe
    df = pd.read_csv(csv, sep='\t')
    df.drop(df.columns[[0]], axis=1,inplace=True)
    #create a plot object
    plt.figure()
    #plot the data
    ax = df.plot(kind='box',xlim=[0,9],ylim=limit,title="{}X Coverage Simulation".format(coverage), grid=True)
    ax.set_xlabel('Loci')
    ax.set_ylabel(label)
    plt.xticks([1,2,3,4,5,6,7,8],['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA'])
    plt.savefig(name)
    plt.close('all')
    
currentPath = os.getcwd()
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--simulationFolder", required = True)
ap.add_argument("-o", "--outputFolderPath", required=True)
args = vars(ap.parse_args())

if not os.path.exists(args["outputFolderPath"]):
    os.mkdir(args["outputFolderPath"])

csvfiles = [f for f in os.listdir("{0}/{1}".format(currentPath, args["simulationFolder"])) if f.endswith(".csv")]

for csv in csvfiles:
    coverage = csv.split("_")[0][:-1]
    if "totalVarDist_bayes" in csv:
        construct_boxPlot(currentPath+"/"+args["simulationFolder"]+ "/"+csv, "TotalVarDist", "{0}/{1}X_totalVarDist_bayes.png".format(args["outputFolderPath"], int(coverage)), int(coverage))
    elif "totalVarDist_count" in csv:
        construct_boxPlot(currentPath+"/"+args["simulationFolder"]+ "/"+csv, "TotalVarDist", "{0}/{1}X_totalVarDist_count.png".format(args["outputFolderPath"], int(coverage)), int(coverage))
    elif "precision" in csv:
        construct_boxPlot(currentPath+"/"+args["simulationFolder"]+ "/"+csv, "Precision", "{0}/{1}X_precision.png".format(args["outputFolderPath"], int(coverage)), int(coverage))
    elif "recall" in csv:
        construct_boxPlot(currentPath+"/"+args["simulationFolder"]+ "/"+csv, "Recall", "{0}/{1}X_recall.png".format(args["outputFolderPath"], int(coverage)), int(coverage))
    elif "diffObjVal" in csv:
        construct_boxPlot(currentPath+"/"+args["simulationFolder"]+ "/"+csv, "DiffObjVal", "{0}/{1}X_diffObjVal.png".format(args["outputFolderPath"], int(coverage)), int(coverage))
        

