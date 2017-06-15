#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 16:14:39 2017

@author: elijah
This scripts similuate a set of strains randomly chosen with random proportions as well.
We are also looking to introduce new strains as errors
"""

from __future__ import division
from collections import defaultdict
from scipy.special import comb
import pandas as pd
import sh
import csv
import math
import numpy as np
import random
import os
import matplotlib.pyplot as plt
import itertools
import sys
import returnStrAndProp as rsp
import argparse

def Generate_Random_strain():
    loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
    random_strain = []
    vals = random.sample(xrange(250,300),len(loci))
    for i in range(len(loci)):
        temp = loci[i]+'_'+str(vals[i])
        random_strain.append(temp)
    return random_strain
        
    
def Compute_Variant_proportions(strain,strain_prop):
    variant_proportions = defaultdict(list)
    for gene in strain:
        variant_proportions[gene] = strain_prop*100
    return variant_proportions
    
def Dict_to_csv(gene_dict, filename):
    with open(filename+'_proportions.csv', 'wb') as csv_file: 
        writer = csv.writer(csv_file)
        for key, value in gene_dict.items():
            writer.writerow([key, value])
            
def Compute_Prec_and_rec(strain_dict, strain_df):
    count = 0
    total_var_dist = 0
    predicted_dict = defaultdict(list)
    for i in range(strain_df.shape[0]):
        row = strain_df.iloc[i].tolist()
        temp_key = tuple(row[0:8])
        predicted_dict[temp_key] = float(row[8])
    print ('The True strains are {}'.format(strain_dict.keys()))
    print ('The True proportions are {}'.format(strain_dict.values()))
    print ('The Predicted strains are {}'.format(strain_dict.keys()))
    print ('The Predicted proportions are {}'.format(strain_dict.values()))
        
    for key in predicted_dict.keys():
        if key in strain_dict.keys():
            
            total_var_dist += abs(round(strain_dict[key],3)-predicted_dict[key])
            count += 1
        else:
            total_var_dist += predicted_dict[key]
    precision = count/strain_df.shape[0]
    recall = count/len(strain_dict)
    return precision, recall, total_var_dist/2.0

def Simulate_strains(numStrains, numIter,strainRef,sampleDir,outputDir,samplesDir):
    
    true_proportions_list = [] #list to hold the true proportions
    true_strains_list = [] #list to hold the true strains
    
    seed = 1994
    random.seed(seed)
    
    
    loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
    precision = []
    recall = []
    total_var_dist = []
    reference = pd.read_csv(strainRef,sep="\t",usecols=range(1,len(loci)+1))
    for name in loci:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
    for iteration in range(1,numIter+1):
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running Simulation {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(iteration))
        #Go to the current sample directory
        os.chdir(sampleDir)
        #define dictionaries to hold each of the variants being used in the simulated strains
        clpA = defaultdict(list)
        clpX = defaultdict(list)
        nifS = defaultdict(list)
        pepX = defaultdict(list)
        pyrG = defaultdict(list)
        recG = defaultdict(list)
        rplB = defaultdict(list)
        uvrA = defaultdict(list)
        #list to hold the names of the varaint dictionaries
        dict_list = [clpA,clpX,nifS,pepX,pyrG,recG,rplB,uvrA]
        #dictionary to hold the strains randomly simulated
        strain_dict = defaultdict(list)
        strains = []
        #generate random proportions for the strains
        proportions = [random.random() for j in range(numStrains)]
        prop_sum = sum(proportions)
        proportions = [i/prop_sum for i in proportions]
        true_proportions_list.append(proportions)
        #pick a random number of unique strains to create
        num_unique = random.randint(1,math.ceil(numStrains/2))
        known = numStrains - num_unique
        #randomly select strains to use for simulation
        randomStrainIndex = random.sample(xrange(1,731),known)
        for num in range(len(randomStrainIndex)):
            index = randomStrainIndex[num]
            strain = reference.iloc[index,:].tolist()
            strains.append(strain)
            tup_strain = tuple(strain)
            strain_dict[tup_strain] = proportions[num]*100
            strain_var_proportions = Compute_Variant_proportions(strain,proportions[num])
            #update the variant dictionaries
            #for clpA
            for key in strain_var_proportions.keys():
                if 'clpA' in key:
                    if key not in clpA:
                        clpA[key] = strain_var_proportions[key]
                    else:
                        clpA[key] = clpA[key] + strain_var_proportions[key]
            #for clpX
            for key in strain_var_proportions.keys():
                if 'clpX' in key:
                    if key not in clpX:
                        clpX[key] = strain_var_proportions[key]
                    else:
                        clpX[key] = clpX[key] + strain_var_proportions[key]
            #for nifS
            for key in strain_var_proportions.keys():
                if 'nifS' in key:
                    if key not in nifS:
                        nifS[key] = strain_var_proportions[key]
                    else:
                        nifS[key] = nifS[key] + strain_var_proportions[key]
            #for pepX
            for key in strain_var_proportions.keys():
                if 'pepX' in key:
                    if key not in pepX:
                        pepX[key] = strain_var_proportions[key]
                    else:
                        pepX[key] = pepX[key] + strain_var_proportions[key]
            #for pyrG
            for key in strain_var_proportions.keys():
                if 'pyrG' in key:
                    if key not in pyrG:
                        pyrG[key] = strain_var_proportions[key]
                    else:
                        pyrG[key] = pyrG[key] + strain_var_proportions[key]
            #for recG
            for key in strain_var_proportions.keys():
                if 'recG' in key:
                    if key not in recG:
                        recG[key] = strain_var_proportions[key]
                    else:
                        recG[key] = recG[key] + strain_var_proportions[key]
            #for rplB
            for key in strain_var_proportions.keys():
                if 'rplB' in key:
                    if key not in rplB:
                        rplB[key] = strain_var_proportions[key]
                    else:
                        rplB[key] = rplB[key] + strain_var_proportions[key]
            #for uvrA
            for key in strain_var_proportions.keys():
                if 'uvrA' in key:
                    if key not in uvrA:
                        uvrA[key] = strain_var_proportions[key]
                    else:
                        uvrA[key] = uvrA[key] + strain_var_proportions[key]
        #generate unique strains and use it for simulation
        for i in range(num_unique):
            strain = Generate_Random_strain()
            strains.append(strain)
            tup_strain = tuple(strain)
            strain_dict[tup_strain] = proportions[i+known]*100
            strain_var_proportions = Compute_Variant_proportions(strain,proportions[i+known])
            for key in strain_var_proportions.keys():
                if 'clpA' in key:
                    if key not in clpA:
                        clpA[key] = strain_var_proportions[key]
                    else:
                        clpA[key] = clpA[key] + strain_var_proportions[key]
            #update the variant dictionaries
            #for clpX
            for key in strain_var_proportions.keys():
                if 'clpX' in key:
                    if key not in clpX:
                        clpX[key] = strain_var_proportions[key]
                    else:
                        clpX[key] = clpX[key] + strain_var_proportions[key]
            #for nifS
            for key in strain_var_proportions.keys():
                if 'nifS' in key:
                    if key not in nifS:
                        nifS[key] = strain_var_proportions[key]
                    else:
                        nifS[key] = nifS[key] + strain_var_proportions[key]
            #for pepX
            for key in strain_var_proportions.keys():
                if 'pepX' in key:
                    if key not in pepX:
                        pepX[key] = strain_var_proportions[key]
                    else:
                        pepX[key] = pepX[key] + strain_var_proportions[key]
            #for pyrG
            for key in strain_var_proportions.keys():
                if 'pyrG' in key:
                    if key not in pyrG:
                        pyrG[key] = strain_var_proportions[key]
                    else:
                        pyrG[key] = pyrG[key] + strain_var_proportions[key]
            #for recG
            for key in strain_var_proportions.keys():
                if 'recG' in key:
                    if key not in recG:
                        recG[key] = strain_var_proportions[key]
                    else:
                        recG[key] = recG[key] + strain_var_proportions[key]
            #for rplB
            for key in strain_var_proportions.keys():
                if 'rplB' in key:
                    if key not in rplB:
                        rplB[key] = strain_var_proportions[key]
                    else:
                        rplB[key] = rplB[key] + strain_var_proportions[key]
            #for uvrA
            for key in strain_var_proportions.keys():
                if 'uvrA' in key:
                    if key not in uvrA:
                        uvrA[key] = strain_var_proportions[key]
                    else:
                        uvrA[key] = uvrA[key] + strain_var_proportions[key]
                
            #strain_var_proportions[strain] = Compute_Variant_proportions(strain,proportions[i+known])
            
        true_strains_list.append(strains)
    
    #write the dictionaries to a csv file
        for i in range(len(dict_list)):
            locus = dict_list[i]
            Dict_to_csv(locus,loci[i])
        
        os.chdir(samplesDir)
        strain_df = rsp.strainSolver(samplesDir,strainRef,outputDir,loci) 
        pre,rec,tvd = Compute_Prec_and_rec(strain_dict,strain_df)
        precision.append(pre)
        recall.append(rec)
        total_var_dist.append(tvd)
        
        
        
        
    
    return precision,recall,total_var_dist
    
    
    


def main():
    #get and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--numOfStrains", required = False, default = 2, type=int, help="Number of strains to simulate. Default = 2")
    ap.add_argument("-n", "--numOfiter", required=False, default=40, type=int, help="The number of simulation iterations default = 40")
    ap.add_argument("-d", "--masterDir", required=False, default=os.getcwd(), help="The master directory to run the simulation and store the results")
    ap.add_argument("-r", "--strainRef", required=True, help="The absolute path to the text file containing the reference strains")
    #ap.add_argument("-p", "--proportionMethod", required=True)
    args = vars(ap.parse_args())
    
    outputDir = "Simulation_Output"
    sampleDir = "Simulation_Samples"
    sample = "Simulation_001"
    
    
    os.chdir(args["masterDir"]) #change into the master directory
    directoriesHere = [d for d in os.listdir(".") if os.path.isdir(d)]
    if sampleDir not in directoriesHere and sample not in directoriesHere and outputDir not in directoriesHere:
        os.mkdir(sampleDir) #make a directory containing all the samples
        os.mkdir(outputDir) #Make a directory for the outputs
        os.chdir(sampleDir)
        os.mkdir(sample) #make a directory for the current simulation sample
    
    #run the simulation
    precision,recall,total_var_dist = Simulate_strains(args["numOfStrains"],args["numOfiter"],args["strainRef"],args["masterDir"]+'/'+sampleDir+"/"+sample,args["masterDir"]+'/'+outputDir,args["masterDir"]+'/'+sampleDir)
    #plot the results
    #for recall
    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creating plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    plt.figure()
    plt.hist(recall, bins=np.linspace(0,2))
    plt.title("Plot of simulation recall")
    plt.xlabel("Recall")
    plt.ylabel("Frequency")
    plt.savefig(args["masterDir"]+'/'+outputDir+'/recall_plot')
    #for precision
    plt.figure()
    plt.hist(precision, bins=np.linspace(0,2))
    plt.title("Plot of simulation precision")
    plt.xlabel("Precision")
    plt.ylabel("Frequency")
    plt.savefig(args["masterDir"]+'/'+outputDir+'/Precision_plot')
    #for total variation distance
    plt.figure()
    plt.hist(total_var_dist,bins=np.linspace(-1,1))
    plt.title("Plot of simulation Total variation distance")
    plt.xlabel("Total Variation Distance")
    plt.ylabel("Frequency")
    plt.savefig(args["masterDir"]+'/'+outputDir+'/Total_variation_Distance_plot')
    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    print ("The recall values for all simulations are: {}".format(recall))
    print ("The precision values for all simulations are: {}".format(precision))
    print ("The total variation distances for all simulations are: {}".format(total_var_dist))
    
#run the main function
if __name__ == "__main__":
    main()
 


        
        
        
        
    
    
    
    
