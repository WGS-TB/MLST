#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 13:11:41 2018

@author: Elijah Willie

Project: Illuminating the diversity of pathogenic bacteria Borrelia Burgdorferi in tick samples

The below program runs the StrainEST algorithm on simulated data based on the new evolution models.
"""

#import the neccessary libraries
from __future__ import division
import pipeline_functions as pf
import pandas as pd
import os
import random
import utility_functions as uf
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
plt.style.use('ggplot')

#set seed for reproducibility 
seed = 1995
random.seed(seed)

#function to run the simulation
def Generate_Data(editDist, num_mut, model, iteration,drop_strain=None):
    
    ''' 
    A function that sets up the type of simuation we are running and generates the neccessary data
    
    Inputs:
        editDist-The maximum edit distance to use
        num_mut-The maximum number of mutations to introduce to a strain
        Model-The type of evolutionary model used to generate the data
        iteration-The current iteration that is bein run
        
    Ouputs:
        A list of the true strains used to generate the data
          
    '''
    
    #get the current directory
    curr_dir = os.getcwd()
    
    #set the loci
    loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
    
    #read in the database of strains
    path = os.path.join(curr_dir,'strain_ref.txt')
    reference = pd.read_csv(path,sep="\t",usecols=range(1,len(loci)+1))
    for name in loci:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
        
    #add all the strains to a list    
    all_strains = []
    for index in range(reference.shape[0]):
        all_strains.append(reference.iloc[index,:].values.tolist())
        
    #if we are running the first evolutionary model
    if model == 'EvoMod_1':
        #run the first evolutionary model
        proportions, true_all, true_alleles, strains = uf.EvoMod1(reference, editDist, num_mut, iteration, all_strains, seed)
        
    #if we are running secondary evolutionary model
    elif model == 'EvoMod_2':
        proportions, true_all, true_alleles, strains = uf.EvoMod2(reference, editDist, num_mut, iteration, all_strains, seed, dropStrain=drop_strain)
        
    #if we are running the third evolutionary model    
    elif model == 'EvoMod_3':
        proportions, true_all, true_alleles, strains = uf.EvoMod3(reference, editDist, num_mut, iteration, all_strains, seed)
    
    #write the strains and their proportions into a dictionary
    strain_prop_dict = uf.Write_Proportions(strains, proportions)
    return strains, strain_prop_dict, true_all
    
def main():
    #get and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--numOfiter", required=False, default=40, type=int, help="The number of simulation iterations default = 40")
    ap.add_argument("-st", "--modType", required=False, default="EvoMod_1", type=str, help="The type of evolution experiment to run. Default = EvoMod_1")
    ap.add_argument("-Ed", "--editDist", required=False, default=15, type=int, help="The maximum edit distance. Default = 15")
    ap.add_argument("-nm", "--numMut", required=False, default=1, type=int, help="The number of mutations to introduce. Default = 1")
    ap.add_argument("-em2v", "--em2v", required=False, default=None, help="'new'/'existing', a variant of evomod2 whether to drop new or existing strain")
    #parse the arguments
    args = vars(ap.parse_args())
    
    print('--------------------------------------------------You have chosen to run {} with {} mutations and edit distance {}-------------------------------------------------------'.format(args['modType'], args['numMut'], args['editDist']))
    
    #get the current directory
    curr_dir = os.getcwd()
    #set the samples directory for the ADP
    samplesDir = os.getcwd()
    #set the Output directory for the ADP
    if not os.path.exists('SDP_Output'):
        os.mkdir('SDP_Output')
    outputDir = os.path.join(curr_dir, 'SDP_Output')
    #set the path the distance matrices for the SDP
    pathToDistMat = os.path.join(curr_dir, 'Dist_DB')
    #set the path to the strain reference 
    path = os.path.join(curr_dir,'strain_ref.txt')
    #set the loci
    loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
    
    SEST_prec = []
    SEST_rec = []
    SEST_softprec = []
    SEST_softrec = []
    SEST_EMD = []
    SEST_TVD = []    

    SEST_Allele_prec = []
    SEST_Allele_rec = []
    SEST_Allele_tvd = []

    ADP_prec_list = []
    ADP_rec_list = []
    ADP_TVD_list = []
    
    SDP_prec_list = []
    SDP_rec_list = []
    SDP_TVD_list = []
    SDP_softprec_list = []
    SDP_softrec_list = []
    SDP_EMD_list = []    

    for i in range(args["numOfiter"]):
        print('-----------------------------------------------------NOW RUNNING SIMULATION {} -----------------------------------------------------------------------------------------'.format(i))
        #generate the data
        strains, strains_prop_dict, true_all = Generate_Data(args["editDist"], args["numMut"], args["modType"], i, args['em2v'])
        
        #compute the alleles
        true_alleles = set([item for sublist in strains for item in sublist])
        
        #run the ADP
        ADP_dict, ADP_alleles, ADP_prop = uf.run_ADP(args["editDist"])
        #compute the ADP statisitcs
        ADP_prec, ADP_rec, ADP_tvd = uf.Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        #append the results
        ADP_prec_list.append(ADP_prec)
        ADP_rec_list.append(ADP_rec)
        ADP_TVD_list.append(ADP_tvd)
        
        
        #run the SDP
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~') 
        strain_df =  pf.strainSolver(samplesDir,path,outputDir,'noProp','s',10800,5, loci=loci,pathToDistMat=pathToDistMat)
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[2:10])
        #compute the SDP statistics
        SDP_prec, SDP_rec, SDP_tvd = uf.Compute_Prec_and_rec(strains_prop_dict, strain_df)
        SDP_softprec, SDP_softrec = uf.Compute_Soft_Prec_and_Rec(strains_prop_dict, strain_df, uf._EDITDIST)
        true_dict = strains_prop_dict
        pred_dict = dict()
        for i in range(strain_df.shape[0]):
            row = strain_df.iloc[i].tolist()
            temp_key = tuple(row[2:10])
            pred_dict[temp_key] = float(row[10])
        #Sort the genes in the strains to be consistent
        true_dict = {tuple(sorted(i)):true_dict[i] for i in true_dict.keys()}
        pred_dict = {tuple(sorted(i)):pred_dict[i] for i in pred_dict.keys()}

        SDP_EMD = uf.compute_EMD(pred_dict, true_dict)
        #append the results
        SDP_prec_list.append(SDP_prec)
        SDP_rec_list.append(SDP_rec)
        SDP_TVD_list.append(SDP_tvd)
        SDP_softprec_list.append(SDP_softprec)
        SDP_softrec_list.append(SDP_softrec)
        SDP_EMD_list.append(SDP_EMD)

        #run StrainEST and get the results
        StrainEST_result, predicted_alleles, result_dict,allele_dict = uf.Run_strainEST(args["editDist"])
        #compute the precision and recall values
        prec, rec, TVD = uf.Compute_SEST_Prec_and_rec(strains_prop_dict, result_dict)
        softprec, softrec = uf.Compute_SEST_Soft_Prec_and_Rec(strains_prop_dict, result_dict, uf._EDITDIST)
        sest_pred_dict = {tuple(sorted(i)):result_dict[i] for i in result_dict.keys()}
        EMD = uf.compute_EMD(true_dict, sest_pred_dict)
        #compute the allele prec and rec
        allele_prec, allele_rec, allele_tvd = uf.Compute_ADP_Prec_and_rec(true_all, allele_dict)
        
        #append these values to their respective lists
        SEST_prec.append(prec)
        SEST_rec.append(rec)
        SEST_TVD.append(TVD)
        SEST_EMD.append(EMD)
        SEST_softprec.append(softprec)
        SEST_softrec.append(softrec)
        SEST_Allele_prec.append(allele_prec)
        SEST_Allele_rec.append(allele_rec)
        SEST_Allele_tvd.append(allele_tvd)
    
    #print the strain results to consle
    print('~~~~~~~~~~~~~~ StrainEST ~~~~~~~~~~~~~~~~~~~~~')
    print('The StrainEST precision values are: {}'.format(SEST_prec))
    print('The StrainEST precision average is: {}'.format(sum(SEST_prec)/len(SEST_prec)))
    print('The StrainEST precision SD is: {}'.format(np.std(np.array(SEST_prec))))
    print('')
    print('The StrainEST Soft precision values are: {}'.format(SEST_softprec))
    print('The StrainEST Soft precision average is: {}'.format(sum(SEST_softprec)/len(SEST_softprec)))
    print('The StrainEST Soft precision SD is: {}'.format(np.std(np.array(SEST_softprec))))
    print('')
    print('The StrainEST recall values are: {}'.format(SEST_rec))
    print('The StrainEST recall average is: {}'.format(sum(SEST_rec)/len(SEST_rec)))
    print('The StrainEST recall SD is: {}'.format(np.std(np.array(SEST_rec))))
    print('')
    print('The StrainEST soft recall values are: {}'.format(SEST_softrec))
    print('The StrainEST soft recall average is: {}'.format(sum(SEST_softrec)/len(SEST_softrec)))
    print('The StrainEST soft recall SD is: {}'.format(np.std(np.array(SEST_softrec))))
    print('')
    print('The StrainEST Total Variation Distances are: {}'.format(SEST_TVD))
    print('The StrainEST Total Variation Distance average is: {}'.format(sum(SEST_TVD)/len(SEST_TVD)))
    print('The StrainEST Total Variation Distance SD is: {}'.format(np.std(np.array(SEST_TVD))))
    print('')
    print('The StrainEST Earth Mover Distances are: {}'.format(SEST_EMD))
    print('The StrainEST Earth Mover Distance average is: {}'.format(sum(SEST_EMD)/len(SEST_EMD)))
    print('The StrainEST Earth Mover Distance SD is: {}'.format(np.std(np.array(SEST_EMD))))
    print('')

     #print allele results to consle
    print('The StrainEST Allele precision is: {}'.format(SEST_Allele_prec))
    print('The StrainEST Allele precision average is: {}'.format(sum(SEST_Allele_prec)/len(SEST_Allele_prec)))
    print('The StrainEST Allele precision SD is: {}'.format(np.std(np.array(SEST_Allele_prec))))
    print('')
    print('The StrainEST Allele recall is: {}'.format(SEST_Allele_rec))
    print('The StrainEST Allele recall average is: {}'.format(sum(SEST_Allele_rec)/len(SEST_Allele_rec)))
    print('The StrainEST Allele recall SD is: {}'.format(np.std(np.array(SEST_Allele_rec))))
    print('')
    print('The StrainEST Allele TVD is: {}'.format(SEST_Allele_tvd))
    print('The StrainEST Allele TVD average is: {}'.format(sum(SEST_Allele_tvd)/len(SEST_Allele_tvd)))
    print('The StrainEST Allele TVD SD is: {}'.format(np.std(np.array(SEST_Allele_tvd))))
    print('')

    print('~~~~~~~~~~~~~~~~~~~~~~~ Our algorithm ~~~~~~~~~~~~~~~~~~~~~')
    #print the ADP results to consle
    print('The ADP precision values are: {}'.format(ADP_prec_list))
    print('The ADP precision average is: {}'.format(sum(ADP_prec_list)/len(ADP_prec_list)))
    print('The ADP precision SD is: {}'.format(np.std(np.array(ADP_prec_list))))
    print('')
    print('The ADP recall values are: {}'.format(ADP_rec_list))
    print('The ADP recall average is: {}'.format(sum(ADP_rec_list)/len(ADP_rec_list)))
    print('The ADP recall SD is: {}'.format(np.std(np.array(ADP_rec_list))))
    print('')
    print('The ADP Total Variation Distances are: {}'.format(ADP_TVD_list))
    print('The ADP Total Variation Distance average is: {}'.format(sum(ADP_TVD_list)/len(ADP_TVD_list)))
    print('The ADP Total Variation Distance SD is: {}'.format(np.std(np.array(ADP_TVD_list))))
    print('')
    
    #print the SDP results to consle
    print('The SDP precision values are: {}'.format(SDP_prec_list))
    print('The SDP precision average is: {}'.format(sum(SDP_prec_list)/len(SDP_prec_list)))
    print('The SDP precision SD is: {}'.format(np.std(np.array(SDP_prec_list))))
    print('')
    print('The SDP Soft precision values are: {}'.format(SDP_softprec_list))
    print('The SDP Soft precision average is: {}'.format(sum(SDP_softprec_list)/len(SDP_softprec_list)))
    print('The SDP Soft precision SD is: {}'.format(np.std(np.array(SDP_softprec_list))))
    print('')
    print('The SDP recall values are: {}'.format(SDP_rec_list))
    print('The SDP recall average is: {}'.format(sum(SDP_rec_list)/len(SDP_rec_list)))
    print('The SDP recall SD is: {}'.format(np.std(np.array(SDP_rec_list))))
    print('')
    print('The SDP Soft recall values are: {}'.format(SDP_softrec_list))
    print('The SDP Soft recall average is: {}'.format(sum(SDP_softrec_list)/len(SDP_softrec_list)))
    print('The SDP Soft recall SD is: {}'.format(np.std(np.array(SDP_softrec_list))))
    print('')
    print('The SDP Total Variation Distances are: {}'.format(SDP_TVD_list))
    print('The SDP Total Variation Distance average is: {}'.format(sum(SDP_TVD_list)/len(SDP_TVD_list)))
    print('The SDP Total Variation Distance SD is: {}'.format(np.std(np.array(SDP_TVD_list))))
    print('')
    print('The SDP Earth Mover Distances are: {}'.format(SDP_EMD_list))
    print('The SDP Earth Mover Distance average is: {}'.format(sum(SDP_EMD_list)/len(SDP_EMD_list)))
    print('The SDP Earth Mover Distance SD is: {}'.format(np.std(np.array(SDP_EMD_list))))
    print('')

    ''' plot the StrainEST statistics '''
    #create a directory to hold the results
    if not os.path.exists(os.path.join(curr_dir, 'Results')):
        os.mkdir(os.path.join(curr_dir, 'Results'))
    #change into that directory
    os.chdir(os.path.join(curr_dir, 'Results'))
    
    #plot the strainEST strain results
    #precision
    plt.figure()
    plt.hist(SEST_prec, bins=np.linspace(0,2))
    plt.title("Plot of StrainEST strain precision for {}".format(args['modType']))
    plt.xlabel("Precision")
    plt.ylabel("Frequency")
    plt.savefig('StrainEST_Strain_Precision_Plot')
    #Save the plot for boxplot plotting
    SEST_precDF = pd.DataFrame(SEST_prec)
    SEST_precDF = SEST_precDF.T
    SEST_precDF.to_csv('StrainEST_Strain_Precision.csv', sep = '\t')
    
    #recall
    plt.figure()
    plt.hist(SEST_rec, bins=np.linspace(0,2))
    plt.title("Plot of StrainEST strain recall for {}".format(args['modType']))
    plt.xlabel("Recall")
    plt.ylabel("Frequency")
    plt.savefig('StrainEST_Strain_Recall_Plot')
    #Save the plot for boxplot plotting
    SEST_recDF = pd.DataFrame(SEST_rec)
    SEST_recDF = SEST_recDF.T
    SEST_recDF.to_csv('StrainEST_Strain_Recall.csv', sep = '\t')
   
    #Soft precision
    SEST_softprecDF = pd.DataFrame(SEST_softprec)
    SEST_softprecDF = SEST_softprecDF.T
    SEST_softprecDF.to_csv('StrainEST_Strain_SoftPrecision.csv', sep = '\t')

    #Soft recall
    SEST_softrecDF = pd.DataFrame(SEST_softrec)
    SEST_softrecDF = SEST_softrecDF.T
    SEST_softrecDF.to_csv('StrainEST_Strain_SoftRecall.csv', sep = '\t')
 
    #EMD
    SEST_emdDF = pd.DataFrame(SEST_EMD)
    SEST_emdDF = SEST_emdDF.T
    SEST_emdDF.to_csv('StrainEST_Strain_EMD.csv', sep = '\t')

    #TVD
    SEST_tvdDF = pd.DataFrame(SEST_TVD)
    SEST_tvdDF = SEST_tvdDF.T
    SEST_tvdDF.to_csv('StrainEST_Strain_TVD.csv', sep = '\t')

    #plot the strainEST alleles results
    #precision
    plt.figure()
    plt.hist(SEST_Allele_prec, bins=np.linspace(0,2))
    plt.title("Plot of StrainEST allele precision for {}".format(args['modType']))
    plt.xlabel("Precision")
    plt.ylabel("Frequency")
    plt.savefig('StrainEST_Allele_Precision_Plot')
    #Save the plot for boxplot plotting
    SEST_Allele_precDF = pd.DataFrame(SEST_Allele_prec)
    SEST_Allele_precDF = SEST_Allele_precDF.T
    SEST_Allele_precDF.to_csv('StrainEST_Allele_Precision.csv', sep = '\t')
    
    #recall
    plt.figure()
    plt.hist(SEST_Allele_rec, bins=np.linspace(0,2))
    plt.title("Plot of StrainEST allele recall for {}".format(args['modType']))
    plt.xlabel("Recall")
    plt.ylabel("Frequency")
    plt.savefig('StrainEST_Allele_Recall_Plot')
    #Save the plot for boxplot plotting
    SEST_Allele_recDF = pd.DataFrame(SEST_Allele_rec)
    SEST_Allele_recDF = SEST_Allele_recDF.T
    SEST_Allele_recDF.to_csv('StrainEST_Allele_Recall.csv', sep = '\t')
    
    #TVD
    SEST_Allele_tvdDF = pd.DataFrame(SEST_Allele_tvd)
    SEST_Allele_tvdDF = SEST_Allele_tvdDF.T
    SEST_Allele_tvdDF.to_csv('StrainEST_Allele_TVD.csv', sep = '\t')
    
    ''' Our algorithm '''
    #precision
    #Save the plot for boxplot plotting
    SDP_precDF = pd.DataFrame(SDP_prec_list)
    SDP_precDF = SDP_precDF.T
    SDP_precDF.to_csv('SDP_Strain_Precision.csv', sep = '\t')
    
    #recall
    #Save the plot for boxplot plotting
    SDP_recDF = pd.DataFrame(SDP_rec_list)
    SDP_recDF = SDP_recDF.T
    SDP_recDF.to_csv('SDP_Strain_Recall.csv', sep = '\t')
   
    #Soft precision
    SDP_softprecDF = pd.DataFrame(SDP_softprec_list)
    SDP_softprecDF = SDP_softprecDF.T
    SDP_softprecDF.to_csv('SDP_Strain_SoftPrecision.csv', sep = '\t')

    #Soft recall
    SDP_softrecDF = pd.DataFrame(SDP_softrec_list)
    SDP_softrecDF = SDP_softrecDF.T
    SDP_softrecDF.to_csv('SDP_Strain_SoftRecall.csv', sep = '\t')
 
    #EMD
    SDP_emdDF = pd.DataFrame(SDP_EMD_list)
    SDP_emdDF = SDP_emdDF.T
    SDP_emdDF.to_csv('SDP_Strain_EMD.csv', sep = '\t')

    #TVD
    SDP_tvdDF = pd.DataFrame(SDP_TVD_list)
    SDP_tvdDF = SDP_tvdDF.T
    SDP_tvdDF.to_csv('SDP_Strain_TVD.csv', sep = '\t')

    #plot the strainEST alleles results
    #precision
    ADP_Allele_precDF = pd.DataFrame(ADP_prec_list)
    ADP_Allele_precDF = ADP_Allele_precDF.T
    ADP_Allele_precDF.to_csv('ADP_Allele_Precision.csv', sep = '\t')
    
    #recall
    ADP_Allele_recDF = pd.DataFrame(ADP_rec_list)
    ADP_Allele_recDF = SEST_Allele_recDF.T
    ADP_Allele_recDF.to_csv('ADP_Allele_Recall.csv', sep = '\t')
    
    #TVD
    ADP_Allele_tvdDF = pd.DataFrame(ADP_TVD_list)
    ADP_Allele_tvdDF = ADP_Allele_tvdDF.T
    ADP_Allele_tvdDF.to_csv('ADP_Allele_TVD.csv', sep = '\t')

    #change back to top directory
    os.chdir('../')
#run the main function
if __name__ == "__main__":
    main()
    
