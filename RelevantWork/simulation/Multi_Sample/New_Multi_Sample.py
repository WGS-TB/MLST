#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 15:53:12 2018

@author: elijah
"""

from __future__ import division
from collections import defaultdict
import pipeline_functions as pf
import pandas as pd
import os
import random
import utility_functions as uf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse
import cs_sim_strainSolver as cs
plt.style.use('ggplot')


#set seed for reproducibility 
seed = 1995
random.seed(seed)


#create a list to hold all the strains

#function to run the simulation

def run_sim(num_samples, num_strains, editDist, num_mut):
    #a dictionary to hold all the strains
    all_strains_dict = defaultdict(dict)

    #create a list to hold the mutated strains
    strains_mut = []
    
    #create a dictionary to hold the precision, recall and TVD values for the ADP for all samples
    ADP_recall = defaultdict(list)
    ADP_precision = defaultdict(list)
    ADP_TVD = defaultdict(list)
    
    #set the loci
    loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
    
    #read in the database of strains
    loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
    path = '/home/parhamgg/Desktop/MLST_CS/MLST/RelevantWork/simulation/Multi_Sample/strain_ref.txt'
    pathToDistMat = '/home/parhamgg/Desktop/MLST_CS/MLST/RelevantWork/simulation/Multi_Sample/Dist_DB'
    reference = pd.read_csv(path,sep="\t",usecols=range(1,len(loci)+1))
    for name in loci:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
        
    #add all the strains to a list    
    all_strains = []
    for index in range(reference.shape[0]):
        all_strains.append(reference.iloc[index,:].values.tolist())
        
    #select 4 random strains from the database
    strains_orig = [random.choice(all_strains) for _ in range(num_strains)]
                    
    #print(strains_orig)
    #print(len(strains_orig))
    
    #create a list to hold all the strains
    total_strains = []
    
    #now for each strain in the orginal, mutate it and add to the total set of strains
    for strain in strains_orig:
        
        #append current strain to the total set
        total_strains.append(strain)
        
        #mutate the current strain
        strain_mut = uf.Mutate_strain(strain, 15, num_mut)
        
        #check to make sure this strain is not already in the database
        flag = True
        while flag:
            if uf.check_strain(strain_mut,all_strains):
                
                #append mutated strain to list
                total_strains.append(strain_mut)
                
                #append the mutated strain to the total set and to the mutated set
                strains_mut.append(strain_mut)
                
                flag = False
            else:
                strain_mut = uf.Mutate_strain(strain, editDist, num_mut)
                flag = True
        
        

    
    #define the samples
    samples = []
    
    #create each sample
    for i in range(1, num_samples+1):
        samples.append('Sample_{}'.format(i))
        
    #change to the samples directory
    os.chdir('/home/parhamgg/Desktop/MLST_CS/MLST/RelevantWork/simulation/Multi_Sample/Samples')
    
    #for each sample in the sample list, create that sample in the samples directory
    for sample in samples:
        if not os.path.exists(sample):
            os.mkdir(sample)
    
    #change back to top directory
    #os.chdir('../')
    
    #create a dictionary to hold all the samples and their strains
    samples_dict = defaultdict(list)
    
    #create a dictionary that holds the strains as keys, and the samples they are in as values
    strains_dict = defaultdict(list)
    
    #for each sample, select two strains at rondom from the total set of strains and assign it to that sample
    #we also assign assign proportions and generate proportions
    #finally we run the ADP on the current sample
    for sample in samples:
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW PROCESSING {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format(sample))
        #enter the currenct sample directory
        os.chdir(sample)
        
        #compute the strains
        strains = [random.choice(total_strains) for _ in range(2)]

        #update the strains dictionary for the current sample
        for strain in strains:
            strains_dict[tuple(strain)].append(sample)
        #assign the strains
        samples_dict[sample].append(strains)
        
        #assign proportions and generate reads
        proportions = uf.Generate_reads(strains, seed, 15)
        
        #Run the ADP on this sample
        ADP_dict, ADP_alleles, ADP_prop = uf.run_ADP(editDist)
        
        #Compute the precision and recall
        #first create a dictionary that contains all the alleles and their proportions from each of the samples
        true_all = defaultdict(int)
        for strain in strains:
            for gene in strain:
                if gene not in true_all:
                    true_all[gene] = proportions[strains.index(strain)]
                else:
                    true_all[gene] = true_all[gene] + proportions[strains.index(strain)]
        #print(true_all)
        
        #generate the strain proportions dictionary
        strain_prop_dict = uf.Write_Proportions(strains, proportions)
        #add it to the dictionary containing all the strains
        all_strains_dict[sample] = strain_prop_dict
        
        #now compute the precsion, recall, and total variation distance
        ADP_pred, ADP_rec, ADP_tvd = uf.Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        
        #update the respective statistical lists
        ADP_precision[sample].append(ADP_pred)
        ADP_recall[sample].append(ADP_rec)
        ADP_TVD[sample].append(ADP_tvd)
        

        
        #change back to top directory
        os.chdir('../')
        
    #compute the average precision, recall and TVD
    avg_prec, avg_rec, avg_tvd = uf.Compute_avg(ADP_precision,ADP_recall, ADP_TVD)
    #print(all_strains_dict)
    print(strains_dict)
    
    #now we run the SDP on all 4 samples
    samplesDir = '/home/parhamgg/Desktop/MLST_CS/MLST/RelevantWork/simulation/Multi_Sample/Samples'
    outputDir = '/home/parhamgg/Desktop/MLST_CS/MLST/RelevantWork/simulation/Multi_Sample/SDP_Results'
    SDP_result = cs.strainSolver(samplesDir,path,outputDir,'noProp','all',10800,5, loci=loci,pathToDistMat=pathToDistMat)
    #print(SDP_result)
    
    #compute the precision, recall, and tvd for the SDP
    S1 = sorted(all_strains_dict.keys())
    S2 = sorted(SDP_result.keys())
    
    #dictionary to hold the precision, recall and tvd values 
    SDP_prec = defaultdict(list)
    SDP_recall = defaultdict(list)
    SDP_tvd = defaultdict(list)
    
    
    #for each sample in the SDP result, compute the summary statistics
    for i in range(len(S1)):
        #print(S1[i])
        #print(S2[i])
        T1 = SDP_result[S2[i]]
        #print(T1)
        #print(all_strains_dict[S1[i]])
        prec, rec, tvd = uf.Compute_Prec_and_rec(all_strains_dict[S1[i]], T1)
        print(prec, rec, tvd)
        SDP_prec[S1[i]].append(prec)
        SDP_recall[S1[i]].append(rec)
        SDP_tvd[S1[i]].append(tvd)
    
    #compute the average summary statistics
    avg_SDP_prec, avg_SDP_rec, avg_SDP_tvd = uf.Compute_avg(SDP_prec,SDP_recall, SDP_tvd)
    
    
    return avg_prec, avg_rec, avg_tvd, avg_SDP_prec, avg_SDP_rec, avg_SDP_tvd
        
    
def main():
    
    #get and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--numOfiter", required=False, default=40, type=int, help="The number of simulation iterations default = 40")
    ap.add_argument("-s", "--numOfstrains", required=False, default=5 ,type=int, help="The number of initial strains to select. Default = 5")
    ap.add_argument("-Ns", "--numOfsample", required=False, default=4, type=int, help="The number of samples to simulate. Default = 4")
    ap.add_argument("-Ed", "--editDist", required=False, default=15, type=int, help="The maximum edit distance. Default = 15")
    ap.add_argument("-nm", "--numMut", required=False, default=1, type=int, help="The number of mutations to introduce. Default = 1")
    #parse the arguments
    args = vars(ap.parse_args())
    
    #containers to hold the ADP statistics
    ADP_avg_prec = []
    ADP_avg_rec = []
    ADP_avg_tvd = []
    
    #containers to hold the SDP statistics
    SDP_avg_prec = []
    SDP_avg_rec = []
    SDP_avg_tvd = []
    
    #set up the simulation iterations
    for i in range(args["numOfiter"]):
        
        print('---------------------------------------------NOW RUNNING ITERATION {}----------------------------------------------------------'.format(i+1))
        #run the simuation
        avg_prec, avg_rec, avg_tvd, avg_SDP_prec, avg_SDP_rec, avg_SDP_tvd = run_sim(args["numOfsample"], args["numOfstrains"], args["editDist"], args["numMut"])
        
        #update the ADP statistics
        ADP_avg_prec.append(avg_prec)
        ADP_avg_rec.append(avg_rec)
        ADP_avg_tvd.append(avg_tvd)
        
        #update the SDP statistics
        SDP_avg_prec.append(avg_SDP_prec)
        SDP_avg_rec.append(avg_SDP_rec)
        SDP_avg_tvd.append(avg_SDP_tvd)
        
    #print some output to the console
    print('-----------------------------------------------------ADP Statistics----------------------------------------------------------------')
    print('Precision is: {}'.format(ADP_avg_prec))
    print('Average precision is: {}'.format(sum(ADP_avg_prec)/len(ADP_avg_prec)))
    print('Recall is: {}'.format(ADP_avg_rec))
    print('Average recall is: {}'.format(sum(ADP_avg_rec)/len(ADP_avg_rec)))
    print('Total Variation Distance is: {}'.format(ADP_avg_tvd))
    print('Average Total Variation Distance is: {}'.format(sum(ADP_avg_tvd)/len(ADP_avg_tvd)))
    
    
    print('-----------------------------------------------------SDP Statistics----------------------------------------------------------------')
    print('Precision is: {}'.format(SDP_avg_prec))
    print('Average precision is: {}'.format(sum(SDP_avg_prec)/len(SDP_avg_prec)))
    print('Recall is: {}'.format(SDP_avg_rec))
    print('Average recall is: {}'.format(sum(SDP_avg_rec)/len(SDP_avg_rec)))
    print('Total Variation Distance is: {}'.format(SDP_avg_tvd))
    print('Average Total Variation Distance is: {}'.format(sum(SDP_avg_tvd)/len(SDP_avg_tvd)))
    
    
    #plot the statisitics
    
    #Create a Results folder
    if not os.path.exists('Results'):
        os.mkdir('Results')
    
    #change into that directory
    os.chdir('Results')
    
    
    #for the ADP
    #ADP Precision
    plt.figure()
    plt.hist(ADP_avg_prec, bins=np.linspace(0,2))
    plt.title("Plot of ADP precision for Multi Sample Simulation")
    plt.xlabel("Precision")
    plt.ylabel("Frequency")
    plt.savefig('ADP_Multi_Sample_Precision_Plot')
    #Save the plot for boxplot plotting
    ADP_precDF = pd.DataFrame(ADP_avg_prec)
    ADP_precDF = ADP_precDF.T
    ADP_precDF.to_csv('ADP_Multi_Sample_Precision.csv', sep = '\t')
    
    #ADP Recall
    plt.figure()
    plt.hist(ADP_avg_rec, bins=np.linspace(0,2))
    plt.title("Plot of ADP Recall for Multi Sample Simulation")
    plt.xlabel("Recall")
    plt.ylabel("Frequency")
    plt.savefig('ADP_Multi_Sample_Recall_Plot')
    #Save the plot for boxplot plotting
    ADP_recDF = pd.DataFrame(ADP_avg_rec)
    ADP_recDF = ADP_recDF.T
    ADP_recDF.to_csv('ADP_Multi_Sample_Recall.csv', sep = '\t')
    
    #ADP TVD
    plt.figure()
    plt.hist(ADP_avg_tvd, bins=np.linspace(-1,1))
    plt.title("Plot of ADP Total Variation Distance for Multi Sample Simulation")
    plt.xlabel("TVD")
    plt.ylabel("Frequency")
    plt.savefig('ADP_Multi_Sample_TVD_Plot')
    #Save the plot for boxplot plotting
    ADP_TVDDF = pd.DataFrame(ADP_avg_tvd)
    ADP_TVDDF = ADP_TVDDF.T
    ADP_TVDDF.to_csv('ADP_Multi_Sample_TVD.csv', sep = '\t')
    
    #for the SDP
    #ADP Precision
    plt.figure()
    plt.hist(SDP_avg_prec, bins=np.linspace(0,2))
    plt.title("Plot of SDP precision for Multi Sample Simulation")
    plt.xlabel("Precision")
    plt.ylabel("Frequency")
    plt.savefig('SDP_Multi_Sample_Precision_Plot')
    #Save the plot for boxplot plotting
    SDP_precDF = pd.DataFrame(SDP_avg_prec)
    SDP_precDF = SDP_precDF.T
    SDP_precDF.to_csv('SDP_Multi_Sample_Precision.csv', sep = '\t')
    
    #ADP Recall
    plt.figure()
    plt.hist(SDP_avg_rec, bins=np.linspace(0,2))
    plt.title("Plot of SDP Recall for Multi Sample Simulation")
    plt.xlabel("Recall")
    plt.ylabel("Frequency")
    plt.savefig('SDP_Multi_Sample_Recall_Plot')
    #Save the plot for boxplot plotting
    SDP_recDF = pd.DataFrame(SDP_avg_rec)
    SDP_recDF = SDP_recDF.T
    SDP_recDF.to_csv('SDP_Multi_Sample_Recall.csv', sep = '\t')
    
    #ADP TVD
    plt.figure()
    plt.hist(SDP_avg_tvd, bins=np.linspace(-1,1))
    plt.title("Plot of SDP Total Variation Distance for Multi Sample Simulation")
    plt.xlabel("TVD")
    plt.ylabel("Frequency")
    plt.savefig('SDP_Multi_Sample_TVD_Plot')
    #Save the plot for boxplot plotting
    SDP_precDF = pd.DataFrame(SDP_avg_tvd)
    SDP_precDF = SDP_precDF.T
    SDP_precDF.to_csv('SDP_Multi_Sample_TVD.csv', sep = '\t')
    
    
    
    


#run the main function
if __name__ == "__main__":
    main()
