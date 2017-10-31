#!/usr/bin/python
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

#Return x with first letter being capital
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

#Return 2D dictionary
def tree():
    return defaultdict(tree)

#Return boolean whether predicted matches true
def predictedCorrectly(predicted, true):
        return set(predicted) == set(true)
    
def compute_QSum(Qmatrix):
    Qsum = Qmatrix.sum()
    Qsum = Qsum.sum()
    return Qsum

def ConvertAscii(Qscore):
    result = []
    for char in Qscore:
        #subtract 33 for predefined offset
        result.append(ord(char)-33)
    return result

def generate_matrix(path):
    var_list = [] #holds the variants
    read_list = [] #holds the reads
    mismatch_list = [] #holds the mismatches
    first = True
    with open(path) as inf:
        for line in inf:
            if first:
                prevParts = line.split('\t')
            else:
                parts = line.split('\t')
            if len(prevParts) > 0 and not first:
                read_list.append(parts[0])
                var_list.append(parts[1])
                mm_1 = re.findall(r'\d+', prevParts[6])
                mm_2 = re.findall(r'\d+',parts[6])
                
                mm_tot = int(mm_1[0]) + int(mm_2[0])
                mismatch_list.append(mm_tot)
            first = not first
            #parts = line.split('\t') # split line into parts
    #        if len(parts) > 1:   # if at least 2 parts/columns
    #            read_list.append(parts[0]) #append reads to a list
    #            var_list.append(parts[1])   #append vars to a list
    #            #mismatch_list.append(parts[2]) #append mismatches to list
    #            MisMatch_pos = re.findall(r'\d+',parts[7])
    #            Converted_scores = ConvertAscii(parts[5])
    #            Qsum = 0
    #            if len(MisMatch_pos) == 1:
    #                Qsum = 0
    #            else:
    #                print MisMatch_pos
    #                for i in range(1,len(MisMatch_pos)):
    #                    sen = int(MisMatch_pos[0])
    #                    pos = int(MisMatch_pos[i])
    #                    tot = sen+pos
    #                    Qsum += Converted_scores[tot-1]
                        
            
    flag = True #makes sure all the previous steps completed successfully
            
            
    if flag is True:
        read_var_dict = defaultdict(list) #dictionary holding all the variants that a read maps to
        read_mismatch_dict = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
        read_index_dict = defaultdict(list) #dictionary holding indices for later use
    
        for i in range(len(read_list)):
            #num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
            num_mismatch = mismatch_list[i]
            #append the appropriate suffix for paired reads
            '''if  i%2 == 0:    
                read_list[i] = read_list[i]+'-1/2'
            else:
                read_list[i] = read_list[i]+'-2/2'
            '''
            
            read_var_dict[read_list[i]].append(var_list[i]) #append all the variants that read read_i maps to
            read_mismatch_dict[read_list[i]].append(num_mismatch) #append all the mismatches that each read_i has when it maps to a variants
            read_index_dict[read_list[i]].append(i) #for testing purposes
    
        var_list = set(var_list) #removes duplicates
        matrix_dict = tree() #creates a 2-D dictionary object later used to generate the read-variant matrix datastructure
     
    #       create a 2D dictionary that contains all the possible combinations of a read with a variant and the number of mismatches.
        for var in var_list:
            for read in read_var_dict:   #key=read name
                temp_var_list = read_var_dict[read]   #list of variants that key maps to 
                if var in temp_var_list:
                    index = temp_var_list.index(var)
                    mismatch = read_mismatch_dict[read][index] #get the number of mismatches
                    #only for perfect matches.
                    if mismatch <= 6:
                    #print val
                        matrix_dict[read][var] = int(mismatch) #add it to the matrix data structure
    
        matrixDF = pd.DataFrame(matrix_dict).T.fillna(-1) #convert 2-D dictionary to a matrix
    return matrixDF

'''
Input: Path to reads.txt file
Output: Matrix with rows=reads, columns=variants and entries=mismatches information
'''
def Compute_Paired_Matrix(path):
    
    first = True
    var_list = [] #holds the variants
    read_list = [] #holds the reads
    mismatch_list = [] #holds the mismatches
    with open(path) as inf:
        
        for line in inf:
            count = 0
            if first:
                prevParts = line.split('\t')
            else:
                parts = line.split('\t')
            if not first and len(prevParts) > 0:
                if parts[3] == '=':
                    read_list.append(parts[0])
                    var_list.append(parts[1])
                    mm_1 = re.findall(r'\d+', prevParts[6])
                    mm_2 = re.findall(r'\d+',parts[6])
                    
                    mm_tot = int(mm_1[0]) + int(mm_2[0])
                    mismatch_list.append(mm_tot)
                else:
                    read = parts[0]+'-{}'.format(str(count))
                    read_list.append(read)
                    var_list.append(parts[1])
                    mm = re.findall(r'\d+',parts[6])
                    mismatch_list.append(mm)
                    count += 1
            first = not first
    flag = True #makes sure all the previous steps completed successfully
                
                
    if flag is True:
        read_var_dict = defaultdict(list) #dictionary holding all the variants that a read maps to
        read_mismatch_dict = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
        read_index_dict = defaultdict(list) #dictionary holding indices for later use
    
        for i in range(len(read_list)):
            #num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
            num_mismatch = mismatch_list[i]
            #append the appropriate suffix for paired reads
            '''if  i%2 == 0:    
                read_list[i] = read_list[i]+'-1/2'
            else:
                read_list[i] = read_list[i]+'-2/2'
            '''
            
            read_var_dict[read_list[i]].append(var_list[i]) #append all the variants that read read_i maps to
            read_mismatch_dict[read_list[i]].append(num_mismatch) #append all the mismatches that each read_i has when it maps to a variants
            read_index_dict[read_list[i]].append(i) #for testing purposes
    
        var_list = set(var_list) #removes duplicates
        matrix_dict = tree() #creates a 2-D dictionary object later used to generate the read-variant matrix datastructure
     
    #       create a 2D dictionary that contains all the possible combinations of a read with a variant and the number of mismatches.
        for var in var_list:
            for read in read_var_dict:   #key=read name
                temp_var_list = read_var_dict[read]   #list of variants that key maps to 
                if var in temp_var_list:
                    index = temp_var_list.index(var)
                    mismatch = read_mismatch_dict[read][index] #get the number of mismatches
                    #only for perfect matches.
                    if mismatch <= 6:
                    #print val
                        matrix_dict[read][var] = int(mismatch) #add it to the matrix data structure
    
        matrixDF = pd.DataFrame(matrix_dict).T.fillna(-1) #convert 2-D dictionary to a matrix
        matrixDF = matrixDF[(matrixDF.T != -1).any()]
        return matrixDF
    
def Compute_Singleton_Matrix(path):
    var_list = [] #holds the variants
    read_list = [] #holds the reads
    mismatch_list = [] #holds the mismatches
    with open(path) as inf:
        count = 0
        for line in inf:
            parts = line.split('\t')
            read = parts[0]# + '-{}'.format(str(count))
            read_list.append(read)
            var_list.append(parts[1])
            mm = re.findall(r'\d+',parts[6])
            mismatch_list.append(int(mm[0]))
            count += 1
    flag = True
    if flag is True:
        read_var_dict = defaultdict(list) #dictionary holding all the variants that a read maps to
        read_mismatch_dict = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
        read_index_dict = defaultdict(list) #dictionary holding indices for later use
    
        for i in range(len(read_list)):
            #num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
            num_mismatch = mismatch_list[i]
            #append the appropriate suffix for paired reads
            '''if  i%2 == 0:    
                read_list[i] = read_list[i]+'-1/2'
            else:
                read_list[i] = read_list[i]+'-2/2'
            '''
            
            read_var_dict[read_list[i]].append(var_list[i]) #append all the variants that read read_i maps to
            read_mismatch_dict[read_list[i]].append(num_mismatch) #append all the mismatches that each read_i has when it maps to a variants
            read_index_dict[read_list[i]].append(i) #for testing purposes
    
        var_list = set(var_list) #removes duplicates
        matrix_dict = tree() #creates a 2-D dictionary object later used to generate the read-variant matrix datastructure
     
    #       create a 2D dictionary that contains all the possible combinations of a read with a variant and the number of mismatches.
        for var in var_list:
            for read in read_var_dict:
                 #key=read name
                temp_var_list = read_var_dict[read]   #list of variants that key maps to 
                if var in temp_var_list:
                    index = temp_var_list.index(var)
                    mismatch = read_mismatch_dict[read][index] #get the number of mismatches
                    #only for perfect matches.
                    if mismatch <= 3:
                    #print val
                        matrix_dict[read][var] = int(mismatch)*2 #add it to the matrix data structure
        matrixDF = pd.DataFrame(matrix_dict).T.fillna(-1) #convert 2-D dictionary to a matrix
        matrixDF = matrixDF[(matrixDF.T != -1).any()]
        return matrixDF

def Generate_Qmatrix(DataPath):
    var_list = []
    read_list = []
    mismatch = []
    Qscores = []
    
    first = True
    with open(DataPath) as inf:
        for line in inf:
            if first:
                prevParts = line.split('\t')
            else:
                parts = line.split('\t')
            if len(prevParts) > 0 and not first:
                read_list.append(parts[0])
                var_list.append(parts[1])
                Converted_scores = ConvertAscii(prevParts[5])
                Converted_scores_2 = ConvertAscii(parts[5])
                MisMatch_pos = re.findall(r'\d+',prevParts[7])
                MisMatch_pos_2 = re.findall(r'\d+',parts[7])
                Qsum = 0
                Qsum2 = 0
                if len(MisMatch_pos) == 0:
                    Qsum = 0
                else:
                    for i in range(1,len(MisMatch_pos)):
                        sen = int(MisMatch_pos[0])
                        pos = int(MisMatch_pos[i])
                        tot = sen+pos
                        Qsum += Converted_scores[tot-1]
                if len(MisMatch_pos_2) == 0:
                    Qsum2 = 0
                else:
                    for i in range(1,len(MisMatch_pos_2)):
                        sen = int(MisMatch_pos_2[0])
                        pos = int(MisMatch_pos_2[i])
                        tot = sen+pos
                        Qsum2 += Converted_scores_2[tot-1]
                Qscores.append(Qsum+Qsum2)
            first = not first
            
    flag = True
    
    if flag is True:
        read_var_dict = defaultdict(list) #dictionary holding all the variants that a read maps to
        read_mismatch_dict = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
        read_index_dict = defaultdict(list) #dictionary holding indices for later use
    
        for i in range(len(read_list)):
            num_mismatch = Qscores[i] #count the number of mismatches for each read
            #append the appropriate suffix for paired reads
            '''if  i%2 == 0:    
                read_list[i] = read_list[i]+'-1/2'
            else:
                read_list[i] = read_list[i]+'-2/2'
            '''
            read_var_dict[read_list[i]].append(var_list[i]) #append all the variants that read read_i maps to
            read_mismatch_dict[read_list[i]].append(num_mismatch) #append all the mismatches that each read_i has when it maps to a variants
            read_index_dict[read_list[i]].append(i) #for testing purposes
    
        var_list = set(var_list) #removes duplicates
        matrix_dict = tree() #creates a 2-D dictionary object later used to generate the read-variant matrix datastructure
     
    #       create a 2D dictionary that contains all the possible combinations of a read with a variant and the number of mismatches.
        for var in var_list:
            for read in read_var_dict:   #key=read name
                temp_var_list = read_var_dict[read]   #list of variants that key maps to 
                if var in temp_var_list:
                    index = temp_var_list.index(var)
                    mismatch = read_mismatch_dict[read][index] #get the number of mismatches
                    #only for perfect matches.
                    #if val == 1:
                    #print val
                    matrix_dict[read][var] = int(mismatch) #add it to the matrix data structure
    
        QmatrixDF = pd.DataFrame(matrix_dict).T.fillna(186) #convert 2-D dictionary to a matrix
        return QmatrixDF

#compute the probability of read of length n mapping to a variant with k mismatches using the binomial distribution/without 
def compute_probability(n, k):
    #NO_BINOM=True means not using binomial coefficient
    if NO_BINOM:
        b=10**10
    else:
        b = comb(n, k, exact=False)
    
    x = math.pow(0.99,(n-k))
    y = math.pow(0.01,k)
    prob = b*x*y
        
    return prob

'''
Input: Dataframe with rows=reads, columns=variants
Output: The proportions of variants (type list) based on counting method
'''
def count_compute_proportions(dataframe):
    prob_list = [0.0]*dataframe.shape[1]
    
#    weightage = dataframe[dataframe.loc[:, dataframe.columns.tolist()] == 0.0].count(axis=0).tolist()
    
    for row in dataframe.itertuples(index=False):
        mmInfo = [i for i in list(row) if i>=0]
        min_mm = min(mmInfo)
        numOfVar_minMm = len([i for i in list(row) if i== min_mm])
        indexOfVar_minMm = [i for i in range(len(list(row))) if list(row)[i] == min_mm]
        
#        weights = [weightage[i] for i in indexOfVar_minMm]
#        weights = [i/sum(weights) for i in weights]
        
        track=0
        for i in indexOfVar_minMm:
#            prob_list[i] += (1/numOfVar_minMm) * weights[track]
            prob_list[i] += (1/numOfVar_minMm)
            track += 1
                
    normalize_term = 1.0/(sum(prob_list))
    prob_list = [100.0*normalize_term * i for i in prob_list]
    return prob_list
        

'''
Input: Dataframe with rows=reads, columns=variants
Output: The proportions of variants (type list)
'''
def bayes_compute_proportions(dataframe):
    #computes the proportion of a set of variants given a set of reads uing probabilistic methods
    prob_list = [] #a list to hold the probabilities
    for row in dataframe.itertuples(index=False):
        temp_list = list(row)
        #compute the probability for each row in the matrix
        for i in range(len(temp_list)):
            if temp_list[i] >= 0:
                temp_list[i] = compute_probability(152,int(temp_list[i]))
            else:
                temp_list[i] = 0
        total = sum(temp_list)
        #solve for k
        #try except just in case when we encounter the weird issue where the decision variable for a predicted variant = 1 but was not output
        try:
            temp_list = [j*(1.0/total) for j in temp_list]
        except ZeroDivisionError:
            print(total)
            print(temp_list)
            
        prob_list.append(temp_list)
    col_sums = [sum(k) for k in zip(*prob_list)]
    total_sum = sum(col_sums)
    prop_list = [100.0*l*(1/total_sum) for l in col_sums]
    return prop_list     

'''
Input: Dataframe with rows=reads, columns=variants
Output: The objective value of these variants, reads which do not map to true variants if any
'''
def compute_true_objective_val(dataframe):
        objective_val = 0
        #To record reads which do not map to true variants
        bad_read = list()
        #Even no limit on mm, there are reads which are covered by predicted variants but not
        #true variants. Identify largest mm , if this case happens, increment by max_mm+1
        max_mm = max(dataframe.max())
        for row in dataframe.itertuples():
            mmInfo_list = [i for i in list(row)[1:] if i >= 0]
            
            #mmInfo_list will be empty if a read does not map back to the true variants
            if TEST_EMPTY_LIST == True:
                if len(mmInfo_list) > 0:
                    objective_val += min(mmInfo_list)   #Increment by the minimum number of mismatches
                else:
                    objective_val+=max_mm + 1
                    bad_read.append(list(row)[0])
                    print(list(row))
            else:
                if len(mmInfo_list) == 0:
                    bad_read.append(list(row)[0])
                    print(list(row))
                    
                objective_val += min(mmInfo_list)
        objective_val += len(dataframe.columns)     #Increment by the number of variants used
        return objective_val, bad_read

#Create a dictionary given keys and values which are lists
def create_dictionary(keys, vals):
        my_dict = dict()
        if len(keys) == len(vals):
            for i in range(len(keys)):
                my_dict[keys[i]] = vals[i]
        return my_dict 

#Compute the difference between predicted and true objective values
def compute_obj_diff(predicted, true):
        diff = predicted - true
        return diff
'''
Input: A dataframe with rows=reads, columns=variants
Output: Negative log likelihood score of this solution
'''    
def compute_likelihood(df):
    numVar = df.shape[1]
    likelihood_list = list()
    max_mm = 6  #because now it is paired
    
    for row in df.itertuples(index=False):
        read = list(row)
        
        temp = list()
        for i in range(numVar):
            if read[i] == -1:   #treat those reads which do not map having mm=max_mm+1
                prob = (0.01)**(max_mm+1) * (0.99)**(152 - max_mm -1)
                temp.append(prob)
            else:
                prob = (0.01)**(read[i]) * (0.99)**(152 - read[i])
                temp.append(prob)
                
        likelihood_list.append( sum(temp) )
    
    #Similar to method in GAML paper
    likelihood_list = [i/(2.0*152*numVar) for i in likelihood_list]
    neg_log_likelihood = [-1.0*np.log10(j) for j in likelihood_list]
    
    score = sum(neg_log_likelihood)
    return score

#def outputDataForML(trueProp, dataMatrix, csvfile):
#    var_mm_probDict = {(var, i):0 for (var,i) in itertools.product(dataMatrix.columns.tolist(), range(7))}
#    variants = dataMatrix.columns.tolist()
#    
#    for row in dataMatrix.itertuples(index=False):
#        temp_index = list()
#        temp_mm = list()
#        for i in range(len(list(row))):
#            if list(row)[i] >= 0:
#                temp_index.append(i)
#                temp_mm.append(list(row)[i])
#                
#        for j in range(len(temp_index)):
#            var_mm_probDict[(variants[temp_index[j]], temp_mm[j])] += 1
#                            
#    for (var,mm) in var_mm_probDict.keys():
#        var_mm_probDict[(var,mm)] = var_mm_probDict[(var,mm)] * compute_probability(152, mm)
#    
#    data = list()
#    for var in variants:
#        temp_list = list()
#        
#        for i in range(7):
#            temp_list.append(var_mm_probDict[(var, i)])
#            
#        temp_list.append(trueProp[var])
#        data.append(temp_list)
#    
#    with open(csvfile, "a") as f:
#        writer = csv.writer(f)
#        writer.writerows(data)

def outputDataForML(dataMatrix, csvfile, varProp_dict):
    mm = range(7)
    matrixForML = list()
    variants = varProp_dict.keys()
    print(variants)
    proportions = varProp_dict.values()
    print(proportions)
    
    track=0
    for v in variants:
        temp_array = list()
        temp = dataMatrix.loc[:, v].value_counts()
        
        for i in mm:
            if i in temp.index:
                temp_array.append(temp[i])
            else:
                temp_array.append(0)
                
        temp_array.append(proportions[track])
        track += 1
        matrixForML.append(temp_array)

    with open(csvfile, "a") as f:
        writer = csv.writer(f)
        writer.writerows(matrixForML)
            
        
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
    numOfOptimalSol = list()
    #Some counts
    predictedCorrect_count = 0
    predCorrect_bool_list = []
    seed = 1994
    random.seed(seed)
    
    #Handling some output files
#    outputFolderPath = "{0}/{1}/".format(originalPath, simulation_result_folder)
    outputResultTxtFile = "{0}/{1}/{2}_output_stats.txt".format(originalPath, simulation_result_folder, gene)
    sys.stdout = open(outputResultTxtFile, "w")        #Write print codes to outputResultTxtFile
    
    #Count the number of variants for this gene
    variantsTxtPath = "{0}/sim_data/{1}/variants.txt".format(originalPath,gene)
    num_variants = sum(1 for line in open(variantsTxtPath))
    
    '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation starts here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    
    for iteration in range(1,numOfIter+1):  
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
            ART_command = "art_illumina -qL 33 -qs 10 -qs2 15 -k 3 -rs {} -q -ss HS25 -sam -i ".format(seed) +file_name+" -p -l 76 -f "+str(covOfThisVar)+" -m 200 -s 10 -o "+simulation_name + ' >/dev/null 2>&1'
            os.system(ART_command)
        
        #Putting the pairs together for all variants
        appendFirst_cmd = "cat *_1.fq > "+str(upperfirst(gene)) + "_"+str(iteration)+"_1.fa" #append all the first of the pairs together
        appendSecond_cmd ="cat *_2.fq > "+str(upperfirst(gene))+"_"+str(iteration)+"_2.fa" #append all the second of the pairs together
        os.system(appendFirst_cmd)
        os.system(appendSecond_cmd)
        ref = upperfirst(gene)+"_bowtie2"
        bowtie2_mapping_cmd = 'bowtie2 -x {0} -a -p 4 -1 ./{1}_{2}_1.fa -2 ./{1}_{2}_2.fa -S {1}_{2}.sam --al {0}_{1}_unpaired.sam '.format(ref, upperfirst(gene), str(iteration))
        os.system(bowtie2_mapping_cmd)
        
        #convert from sam to bam file
        convert_cmd = 'samtools view -h -b -S {0}_{1}.sam > {0}_{1}.bam'.format(upperfirst(gene),str(iteration))
        filter_cmd = 'samtools view -b -F 4 {0}_{1}.bam > {0}_{1}_mapped.bam'.format(upperfirst(gene),str(iteration))
        paired_cmd = 'samtools view -b -F8 {0}_{1}_mapped.bam > {0}_{1}_paired.bam'.format(upperfirst(gene),str(iteration))
        singleton_cmd = 'samtools view -b -f8 {0}_{1}_mapped.bam > {0}_{1}_singleton.bam'.format(upperfirst(gene),str(iteration))
        paired_parse_cmd = '''samtools view {0}_{1}_paired.bam | awk '{{print $1 "\\t" $3 "\\t" $4 "\\t" $7 "\\t" $8 "\\t" $11 "\\t" $15 "\\t" $19}}' > {0}_{1}_paired_reads.txt'''.format(upperfirst(gene),str(iteration))
        singleton_parse_cmd = '''samtools view {0}_{1}_singleton.bam | awk '{{print $1 "\\t" $3 "\\t" $4 "\\t" $7 "\\t" $8 "\\t" $11 "\\t" $14 "\\t" $18}}' > {0}_{1}_singleton_reads.txt'''.format(upperfirst(gene),str(iteration))
        #run the commands
        os.system(convert_cmd)
        os.system(filter_cmd)
        os.system(paired_cmd)
        os.system(singleton_cmd)
        os.system(paired_parse_cmd)
        os.system(singleton_parse_cmd)
        
        #Remove unneccessary files for the next iteration.    
        os.system("rm {}*".format(gene))
        os.system("rm *.bam")
        os.system("rm *.sam")
        #os.system("rm {}*.sam".format(upperfirst(gene)))
        #os.system("rm {}*.bam".format(upperfirst(gene)))
        #Keep track of true variants
        true_ratios_list.append(fractions)
        true_variants_list.append(true_variants)
        for j in range(0,len(true_variants)):
            key = true_variants[j]
            true_prop[key] = float(fractions[j])*100
        
        #Create data matrix where rows=reads and columns=variants
        paired_readsTxt_path = upperfirst(gene)+ '_'+str(iteration)+'_paired_reads.txt'
        singleton_readsTxt_path = upperfirst(gene)+ '_'+str(iteration)+'_singleton_reads.txt'
        pairedDF = Compute_Paired_Matrix(paired_readsTxt_path)
        pairedDF = generate_matrix(paired_readsTxt_path)
        singletonDF = Compute_Singleton_Matrix(singleton_readsTxt_path)
        dataMatrix = pd.concat([pairedDF,singletonDF])
        dataMatrix = dataMatrix.fillna(-1)
        dataMatrix.rename(columns={'Unnamed: 0': 'Read'}, inplace=True)

        #Qmatrix = Generate_Qmatrix(readsTxt_path)
        #Run the ILP solver
        pred_object_val,var_predicted,reads_cov,all_solutions, all_objective = varSolver.solver(dataMatrix)
        
        
        '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Statistics and Calculations start here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
        
        print "======================================== SIMULATION " + str(iteration) + " ====================================================" +  "\n"    
        
        #Keep the likelihood scores of all optimal solutions
#        minVar_solutions = [sol for sol in all_solutions if len(sol) == min(map(len,all_solutions))]
#        all_solutions = minVar_solutions
        score_list = list()
        min_score = sys.maxint
        print("Number of optimal solutions: {}".format(len(all_solutions)))
        numOfOptimalSol.append(len(all_solutions))
        print all_solutions
        Qscore_list = [] #list to hold all the average Qsocres for the solutions
        
        #Compute negative log likelihood score for each solution
        for i in range(len(all_solutions)):
#            print("Solution:{}".format(all_solutions[i]))
#            print("Objective value: {}".format(all_objective[i]))
#            print("Proportions:{}".format(compute_proportions(dataMatrix.loc[reads_cov, all_solutions[i]])))
            score = compute_likelihood(dataMatrix.loc[reads_cov, all_solutions[i]])
            score_list.append(score)
            
            
            if score <= min_score:
                min_score = score
                
            Qscore_list.append(compute_QAvg(Qmatrix.loc[reads_cov,all_solutions[i]]))
        min_Qscore = np.argmin(Qscore_list)
        var_predicted = all_solutions[min_Qscore]
            
        #If there is one solution among all optimal which matches true variants, assign to var_predicted to generate statistics
        '''for solution in all_solutions:
            if set(solution) == set(true_variants):
                var_predicted = solution
                break
        '''
            
        #Keep track of precision and recall for all simulations
        precision_list.append(precision(var_predicted, true_variants))
        recall_list.append(recall(var_predicted, true_variants))
        pred_object_vals.append(pred_object_val)
        
        #Construct dataframe of true variants and predicted variants
        true_DF = dataMatrix.loc[reads_cov,true_variants]
#        true_DF =  true_DF[(true_DF.T != -1).any()]
        predicted_DF = dataMatrix.loc[reads_cov,var_predicted]
        prop_count = count_compute_proportions(predicted_DF)
        #prop_bayes = bayes_compute_proportions(predicted_DF)
#        if proportion_method == "count":
#            prop = count_compute_proportions(predicted_DF)
#        elif proportion_method == "bayes":
#            prop = bayes_compute_proportions(predicted_DF)
        pred_prop_count = create_dictionary(var_predicted, prop_count)
        #pred_prop_bayes = create_dictionary(var_predicted, prop_bayes)
        val_count = totalVariationDist(pred_prop_count, true_prop)
        #val_bayes = totalVariationDist(pred_prop_bayes, true_prop)
        totalVarDist_count.append(val_count)
        #totalVarDist_bayes.append(val_bayes)
        
        #Print true and predicted variants
        print("True variants are: {}\n".format(true_variants))
        print("Predicted variants are: {}\n".format(var_predicted))
        print("True proportions are: {}\n".format(true_prop))
        #print("Predicted proportions using Bayes method are: {}\n".format(pred_prop_bayes))
        print("Predicted proportions using Counting method are: {}\n".format(pred_prop_count))
        print ('average Qscores are: {}').format(Qscore_list)
        
        #Compute objective value of true variants
        true_Objective_val, bad_reads = compute_true_objective_val(true_DF)
        true_Objective_vals.append(true_Objective_val)
        
        #Compute the difference in objective vlaue: Predicted - True
        diff_obj_vals.append(compute_obj_diff(pred_object_val,true_Objective_val))
        
        #Count simulations in which the ILP predicted correctly
        if predictedCorrectly(var_predicted, true_variants):
                predictedCorrect_count += 1
                predCorrect_bool_list.append(True)
        else:
                predCorrect_bool_list.append(False)
                
        if len(bad_reads) != 0:
            print(dataMatrix.loc[bad_reads,:])   

    
    '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Here is the summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    
    print "======================================== {0}: SUMMARY STATISTICS ====================================================\n".format(gene)
#    avg_totalVarDist_bayes = sum(totalVarDist_bayes)/len(totalVarDist_bayes)
#    true_avg_totalVarDist_bayes = sum(list(itertools.compress(totalVarDist_bayes, predCorrect_bool_list)))/sum(predCorrect_bool_list)
#    variance_totalVarDist_bayes = map(lambda x: (x - avg_totalVarDist_bayes)**2, totalVarDist_bayes)
#    true_variance_totalVarDist_bayes = map(lambda x:(x - true_avg_totalVarDist_bayes)**2, list(itertools.compress(totalVarDist_bayes, predCorrect_bool_list)))
#    variance_totalVarDist_bayes = sum(variance_totalVarDist_bayes)/len(variance_totalVarDist_bayes)
#    true_variance_totalVarDist_bayes = sum(true_variance_totalVarDist_bayes)/len(true_variance_totalVarDist_bayes)
#    std_totalVarDist_bayes = math.sqrt(variance_totalVarDist_bayes)
#    true_std_totalVarDist_bayes = math.sqrt(true_variance_totalVarDist_bayes)
    
    avg_totalVarDist_count = sum(totalVarDist_count)/len(totalVarDist_count)
    true_avg_totalVarDist_count = sum(list(itertools.compress(totalVarDist_count, predCorrect_bool_list)))/sum(predCorrect_bool_list)
    variance_totalVarDist_count = map(lambda x: (x - avg_totalVarDist_count)**2, totalVarDist_count)
    true_variance_totalVarDist_count = map(lambda x:(x - true_avg_totalVarDist_count)**2, list(itertools.compress(totalVarDist_count, predCorrect_bool_list)))
    variance_totalVarDist_count = sum(variance_totalVarDist_count)/len(variance_totalVarDist_count)
    true_variance_totalVarDist_count = sum(true_variance_totalVarDist_count)/len(true_variance_totalVarDist_count)
    std_totalVarDist_count = math.sqrt(variance_totalVarDist_count)
    true_std_totalVarDist_count = math.sqrt(true_variance_totalVarDist_count)
    
    context=1
    print "({0})Total Variation Distance:\n".format(context)
    print "Counting Method ~ Total variation distances are:",totalVarDist_count, "\n"
    #print "Bayes' Method ~ Total variation distances are:",totalVarDist_bayes, "\n"
    print "Counting Method ~ The average of total variation distance is:", avg_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The average of total variation distance is:", avg_totalVarDist_bayes, "\n"
    #print "The variance of total variation distance is:", variance_totalVarDist, "\n"
    print "Counting Method ~ The standard deviation of total variation distance is:",std_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The standard deviation of total variation distance is:",std_totalVarDist_bayes, "\n"
    context+=1
    
    print "({0})Total Variation Distance for variants which are predicted correctly:\n".format(context)
    print "Counting Method ~ Total variation distances are:",list(itertools.compress(totalVarDist_count, predCorrect_bool_list)), "\n"
    #print "Bayes' Method ~ Total variation distances are:",list(itertools.compress(totalVarDist_bayes, predCorrect_bool_list)), "\n"
    print "Counting Method ~ The average of total variation distance is:", true_avg_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The average of total variation distance is:", true_avg_totalVarDist_bayes, "\n"
    #print "The variance of total variation distance is:", variance_totalVarDist, "\n"
    print "Counting Method ~ The standard deviation of total variation distance is:",true_std_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The standard deviation of total variation distance is:",true_std_totalVarDist_bayes, "\n"
    context+=1
    
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
    
    avg_diffObjVal = sum(diff_obj_vals)/len(diff_obj_vals)
    std_diffObjVal = np.std(np.array(diff_obj_vals))
    print "({0}) Objective Value: \n".format(context)
    print 'Predicted objective values are:', pred_object_vals, "\n"
    print 'True objective values are:', true_Objective_vals, "\n"
    print 'The difference in objective values are:', diff_obj_vals, "\n"
    print "Average of difference in objective value is: ", avg_diffObjVal, "\n"
    print "Standard deviation of difference in objective value is: ", std_diffObjVal, "\n"
    context+=1
    
    print "({0})Accuracy: \n".format(context)
    print 'Total number of simulations: ', numOfIter , "\n"
    print("Number of simulations which are predicted correctly: {0}\n".format(predictedCorrect_count))
    print 'Percentage of simulations predicted correctly: ', 100*predictedCorrect_count/iteration, "%\n"
    context+=1
    
    print "({0})Optimal solutions: \n".format(context)
    print("Number of optimal solutions: {0}\n".format(numOfOptimalSol))
    print("Average number of optimal solutions: {0}\n".format(np.mean(numOfOptimalSol)))
    context+=1
    
    sys.stdout.close()
    return precision_list, recall_list, diff_obj_vals, totalVarDist_count