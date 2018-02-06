#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:02:51 2017

@author: elijah
"""

from __future__ import division
from collections import defaultdict
from scipy.misc import comb
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
import variantILP as varSolver
import linecache
import re

#Testing purposes and global variables
NO_BINOM = False
TEST_EMPTY_LIST = True

''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Function Definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

'''
Input: pred_prop and true_prop which are dictionaries
Return: Total variation distance
'''
def totalVariationDist(pred_prop, true_prop):
    common_variants = set.intersection(set(pred_prop.keys()), set(true_prop.keys()))
    diff_vars_pred = set.difference(set(pred_prop.keys()), common_variants)
    diff_vars_true = set.difference(set(true_prop.keys()), common_variants)
    common_variants = list(common_variants)
    diff_vars_pred = list(diff_vars_pred)
    diff_vars_true = list(diff_vars_true)
    
    totalVarDist=0
    for var in common_variants:
        totalVarDist += abs(pred_prop[var] - true_prop[var])
        
    for var in diff_vars_pred:
        totalVarDist += pred_prop[var]
        
    for var in diff_vars_true:
        totalVarDist += true_prop[var]
        
    return totalVarDist/2

def upperfirst(x):
    return x[0].upper() + x[1:]

#predicted and true are lists
def precision(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(predicted))
       
#predicted and true are lists 
def recall(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(true))
#Return boolean whether predicted matches true
def predictedCorrectly(predicted, true):
        return set(predicted) == set(true)

def create_dictionary(keys, vals):
    my_dict = dict()
    if len(keys) == len(vals):
        for i in range(len(keys)):
            my_dict[keys[i]] = vals[i]
    return my_dict 

def simulation(gene, numOfIter, originalPath, simulation_result_folder, coverage):
    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Defining some parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '''
    #Record true variants and their fractions for all simulations
    true_ratios_list = []
    true_variants_list = []
    #Statistics list to output
    totalVarDist_count = []
    totalVarDist_bayes = []
    precision_list = []
    recall_list = []
    pred_object_vals = []
    true_Objective_vals = []
    diff_obj_vals = []
    likelihoodCalibration = []
    numOfOptimalSol = list()
    #Some counts
    predictedCorrect_count = 0
    predCorrect_bool_list = []
#    minNegLogLike_correct = 0
#    minSizeOpt_count = 0
#    minNegLogLike_hasTrue = 0
    #Reproducible results, set seed
    seed = 1994
    random.seed(seed)
    
    #Handling some output files
#    outputFolderPath = "{0}/{1}/".format(originalPath, simulation_result_folder)
    outputResultTxtFile = "{0}/{1}/{2}_output_stats.txt".format(originalPath, simulation_result_folder, gene)
    sys.stdout = open(outputResultTxtFile, "w")        #Write print codes to outputResultTxtFile
    
    #Count the number of variants for this gene
    variantsTxtPath = "{0}/sim_data/{1}/variants.txt".format(originalPath,gene)
    num_variants = sum(1 for line in open(variantsTxtPath))
    
    #Randomly choose some simulations to plot the likelihood graphs
#    randomNum_negLogCharts = random.sample(xrange(1, numOfIter), 10)
    
    '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation starts here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    ref = upperfirst(gene)+'.idx'
    
    for iteration in range(1,numOfIter+1):  
        print "======================================== SIMULATION " + str(iteration) + " ====================================================" +  "\n"    
        true_prop = dict()#dictionary to store the true proportions
        k =random.randint(2,7) #generate a random integer k between 2 and 7
        #generate k random fractions that sum up to 1
        fractions = [random.random() for j in range(k)] 
        s = sum(fractions)
        fractions = [ i/s for i in fractions ]
        #print fractions
        true_variants = []
        #randomly select variants to use for simulated data
        randomVarIndex = random.sample(xrange(1,num_variants+1), k) #start from 1 to num_variants+1 because 0 index of linecache is ''
        for index in randomVarIndex:
            variant =linecache.getline(variantsTxtPath,index) #extract the random variant
            #print variant
            variant = str(variant) 
            variant = variant.rstrip() #remove the "\n" character that is returned by bash
            string1= variant.split(">")
            true_variants.append(string1[1]) #append the variant to the list          
        #print num
        #print true_variants
        
        '''======== Generate the reads using ART ========'''
        total_variants_sequence = ''
        for i in range(len(fractions)):
            variant = true_variants[i]
            simulation_name = true_variants[i] + '_' + str(iteration)+'_'
            file_name = simulation_name + '_reference.fa'
            covOfThisVar = math.ceil((fractions[i]*coverage)) #compute the number of reads to generate
            covOfThisVar = int(covOfThisVar)
            variant_sequence = sh.grep(variant,"{0}/sim_data/{1}/linear.txt".format(originalPath, gene),"-w","-A1") #use bash to extract the variant sequence
            variant_sequence = variant_sequence.rstrip() #remove the "\n" character that is returned by bash
            variant_sequence = str(variant_sequence)
            total_variants_sequence += variant_sequence
            
            #Write sequence to file
            with open(file_name, "w") as sequence_file: 
                sequence_file.write(variant_sequence)
            
            #Set the ART command, I have included a random seed for reproducibility, and a coverage parameter
            ART_command = "art_illumina -qL 33 -qs -10 -k 3 -rs {} -q -ss HS25 -sam -i ".format(seed) +file_name+" -p -l 76 -f "+str(covOfThisVar)+" -m 200 -s 10 -o "+simulation_name + ' >/dev/null 2>&1'
            os.system(ART_command)
        
        #Putting the pairs together for all variants
        appendFirst_cmd = "cat *_1.fq > "+str(upperfirst(gene)) + "_"+str(iteration)+"_1.fa" #append all the first of the pairs together
        appendSecond_cmd ="cat *_2.fq > "+str(upperfirst(gene))+"_"+str(iteration)+"_2.fa" #append all the second of the pairs together
        os.system(appendFirst_cmd)
        os.system(appendSecond_cmd)
        #run kallisto command
        os.system("rm {}*".format(gene))
        
        Kallisto_cmd = 'kallisto quant -t 8 -i {0} -o simulation_{1} ./{2}_{1}_1.fa ./{2}_{1}_2.fa >/dev/null 2>&1'.format(ref,str(iteration),upperfirst(gene))
        os.system(Kallisto_cmd)
        os.chdir('simulation_{}'.format(str(iteration)))
        output_file = pd.read_csv('abundance.tsv',sep='\t')
        DF = output_file.loc[:,['target_id','est_counts']]
        DF = DF[DF['est_counts'] != 0]
        DF['est_counts'] = (DF['est_counts']/DF['est_counts'].sum())*100
        print DF
        DF = DF[DF['est_counts'] > 1.0]
        DF['est_counts'] = (DF['est_counts']/DF['est_counts'].sum())*100
        print "Next"
        print DF
        var_predicted = DF['target_id'].tolist()
        props = DF['est_counts'].tolist()
        #Remove unneccessary files for the next iteration.    
        
        
        #Keep track of true variants
        true_ratios_list.append(fractions)
        true_variants_list.append(true_variants)
        for j in range(0,len(true_variants)):
            key = true_variants[j]
            true_prop[key] = float(fractions[j])*100
        
        pred_dict = create_dictionary(var_predicted,props)
        tvd = totalVariationDist(pred_dict,true_prop)
        totalVarDist_count.append(tvd)
        precision_list.append(precision(var_predicted, true_variants))
        recall_list.append(recall(var_predicted, true_variants))
        
        
        #Print true and predicted variants
        print("True variants are: {}\n".format(true_variants))
        print("Predicted variants are: {}\n".format(var_predicted))
        print("True proportions are: {}\n".format(true_prop))
        print("Predicted proportions are: {}\n".format(pred_dict))
        
        if predictedCorrectly(var_predicted, true_variants):
                predictedCorrect_count += 1
                predCorrect_bool_list.append(True)
        else:
                predCorrect_bool_list.append(False)
        
        #go up a directory
        os.chdir('..')
        
        
        '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Here is the summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    
    print "======================================== {0}: SUMMARY STATISTICS ====================================================\n".format(gene)
   
    
    avg_totalVarDist_count = sum(totalVarDist_count)/len(totalVarDist_count)
    #true_avg_totalVarDist_count = sum(list(itertools.compress(totalVarDist_count, predCorrect_bool_list)))/sum(predCorrect_bool_list)
    variance_totalVarDist_count = map(lambda x: (x - avg_totalVarDist_count)**2, totalVarDist_count)
    #true_variance_totalVarDist_count = map(lambda x:(x - true_avg_totalVarDist_count)**2, list(itertools.compress(totalVarDist_count, predCorrect_bool_list)))
    variance_totalVarDist_count = sum(variance_totalVarDist_count)/len(variance_totalVarDist_count)
    #true_variance_totalVarDist_count = sum(true_variance_totalVarDist_count)/len(true_variance_totalVarDist_count)
    std_totalVarDist_count = math.sqrt(variance_totalVarDist_count)
    #true_std_totalVarDist_count = math.sqrt(true_variance_totalVarDist_count)
    
    context=1
    print "({0})Total Variation Distance:\n".format(context)
    print "Counting Method ~ Total variation distances are:",totalVarDist_count, "\n"
    print "Counting Method ~ The average of total variation distance is:", avg_totalVarDist_count, "\n"
    #print "The variance of total variation distance is:", variance_totalVarDist, "\n"
    print "Counting Method ~ The standard deviation of total variation distance is:",std_totalVarDist_count, "\n"
    context+=1
    
   # print "({0})Total Variation Distance for variants which are predicted correctly:\n".format(context)
   # print "Counting Method ~ Total variation distances are:",list(itertools.compress(totalVarDist_count, predCorrect_bool_list)), "\n"
    #print "Counting Method ~ The average of total variation distance is:", true_avg_totalVarDist_count, "\n"
    #print "The variance of total variation distance is:", variance_totalVarDist, "\n"
   # print "Counting Method ~ The standard deviation of total variation distance is:",true_std_totalVarDist_count, "\n"
    #context+=1
    
    avg_prec = sum(precision_list)/len(precision_list)
    std_prec = np.std(np.array(precision_list))
    print "({0}) Precision: \n".format(context)
    print 'Precision is:', precision_list, "\n"
    print "Average of precision is: ", avg_prec, "\n"
    print "Standard deviation of precision is: ", std_prec, "\n"
    context+=1
    
    avg_rec = sum(recall_list)/len(recall_list)
    std_rec = np.std(np.array(recall_list))
    print "({0}) Recall : \n".format(context)
    print 'Recall is:', recall_list, "\n"
    print "Average of recall is: ", avg_rec, "\n"
    print "Standard deviation of recall is: ", std_rec, "\n"
    context+=1
    
    
    print "({0})Accuracy: \n".format(context)
    print 'Total number of simulations: ', numOfIter , "\n"
    print("Number of simulations which are predicted correctly: {0}\n".format(predictedCorrect_count))
    print 'Percentage of simulations predicted correctly: ', 100*predictedCorrect_count/iteration, "%\n"
    context+=1
    
    sys.stdout.close()
    return precision_list, recall_list, totalVarDist_count
