#!/usr/bin/python
'''
This script take reads(text file) specific to a sample as input, gene name.
It outputs identified variants at that gene and their proportions.
'''
from __future__ import division
from collections import defaultdict
from scipy.special import comb
import pandas as pd
import math
import argparse
import csv
import variantILP as varSolver
import matplotlib.pyplot as plt
import numpy as np
import sys

NO_BINOM = False
'''
Input: Path to reads.txt file
Output: Matrix with rows=reads, columns=variants and entries=mismatches information
'''
def generate_matrix(path):
    var_list = [] #holds the variants
    read_list = [] #holds the reads
    mismatch_list = [] #holds the mismatches
    first=True      #first read
    all_var = set()     #the set of all variants
    
    #Split columns and append them into lists
    with open(path) as inf:
        for line in inf:
            if first:
                prevParts = line.split('\t')
            else:
                parts = line.split('\t')
            if len(prevParts)>1 and not first:
                temp = parts[0].split('-')
                all_var.update([temp[0]])
                temp2 = '-'.join(parts[0].split('-',2)[:2])
                temp2 = temp2.split('/')
                read_list.append(temp2[0])
                var_list.append(parts[1])
                mismatch_list.append(prevParts[2] + parts[2])
            first = not first
        flag = True #makes sure all the previous steps completed successfully
        
        
    if flag is True:
        read_var_dict = defaultdict(list) #dictionary holding all the variants that a read maps to
        read_mismatch_dict = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
        read_index_dict = defaultdict(list) #dictionary holding indices for later use

        for i in range(len(read_list)):
            num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
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

        matrixDF = pd.DataFrame(matrix_dict).T.fillna(-1) #convert 2-D dictionary to a matrix
    return matrixDF

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

#Return 2D dictionary
def tree():
    return defaultdict(tree)

'''
Input: Dataframe with rows=reads, columns=variants
Output: The proportions of variants (type list)
'''
def compute_proportions(dataframe):
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

#Create a dictionary given keys and values which are lists
def create_dictionary(keys, vals):
    my_dict = dict()
    if len(keys) == len(vals):
        for i in range(len(keys)):
            my_dict[keys[i]] = vals[i]
    return my_dict 

'''
Input: A dataframe with rows=reads, columns=variants
Output: Negative log likelihood score of this solution
'''    
def compute_likelihood(df):
    numVar = df.shape[1]
    likelihood_list = list()
    max_mm = 3
    
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

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--path", required = True, help = "Path to data table")
ap.add_argument("-g", "--gene", required = True,  help="Name of gene")
args = vars(ap.parse_args())

path = args["path"]

#generate matrix
dataMatrixDF = generate_matrix(path)
dataMatrixDF.rename(columns={'Unnamed: 0': 'Read'}, inplace=True)
#predict variants
pred_object_val,var_predicted,reads_cov, all_solutions, all_objective = varSolver.solver(dataMatrixDF)
#print var_predicted
#print all_solutions
dataMatrix_pred = dataMatrixDF.loc[reads_cov,var_predicted]

#write matrix to file
dataMatrix_pred.to_csv(args["gene"]+'_predicted_matrix.csv', sep='\t')
fig, ax = plt.subplots()
dataMatrix_pred.hist(bins=np.histogram(dataMatrix_pred.values.ravel())[1],ax=ax)
fig.savefig(args["gene"]+'_matrix_plots.png')

#compute proportions
prop = compute_proportions(dataMatrix_pred)
pred_prop = create_dictionary(var_predicted, prop)

#score list and proportions
score_list = list()
min_score = sys.maxint
min_sol_list = list()
#print("Total number of optimal solutions: {}".format(len(all_solutions)))

for i in range(len(all_solutions)):
#    print("Solution: {}".format(all_solutions[i]))
#    print("Objective value: {}".format(all_objective[i]))
#    print("Proportion: {}".format(compute_proportions(df.loc[reads_cov, all_solutions[i]])))
    score = compute_likelihood(dataMatrixDF.loc[reads_cov, all_solutions[i]])
    score_list.append(score)
    
    if score <= min_score:
        min_score = score
    
#    print("Negative log likelihood score:{}\n".format(score))
    
min_sol_list = [all_solutions[i] for i in range(len(all_solutions)) if score_list[i] == min_score]
    
#print("**Summary**")
#print("Negative log likelihood score list:{}".format(score_list))
#print("Minimum negative log likelihood score: {}".format(min_score))
#print("Solution(s) which have minimum negative log likelihood: {}".format(min_sol_list))

#write proportions to file
w = csv.writer(open(args["gene"]+'_proportions.csv', "w"))
for key, val in pred_prop.items():
    w.writerow([key, val])