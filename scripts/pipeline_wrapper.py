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
import variantILP as ilp
import matplotlib.pyplot as plt
import numpy as np
import sys

#Generate matrix of rows=reads and columns=number of mismatches
def Generate_Matrix(path):
    var_list = [] #holds the variants
    read_list = [] #holds the reads
    mismatch_list = [] #holds the mismatches
    with open(path) as inf:
            for line in inf:
                    parts = line.split('\t') # split line into parts
                    if len(parts) > 1:   # if at least 2 parts/columns
                            read_list.append(parts[0]) #append reads to a list
                            var_list.append(parts[1])   #append vars to a list
                            mismatch_list.append(parts[2]) #append mismatches to list
            flag = True
    if flag is True:
            d = defaultdict(list) #dictionary holding all the variants that a read maps to
            d_1 = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
            d_2 = defaultdict(list) #dictionary holding indices for later use


            for i in range(len(read_list)):
                    num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
                    
                    if  i%2 == 0:    
                        read_list[i] = read_list[i]+'-1/2'
                    else:
                        read_list[i] = read_list[i]+'-2/2'
                    
                    d[read_list[i]].append(var_list[i]) #append all the variants that read read_i maps to
                    d_1[read_list[i]].append(num_mismatch) #append all the mismatches that each read_i has when it maps to a variants
                    d_2[read_list[i]].append(i) #for testing purposes

            var_list = set(var_list)
            x = tree()

#               create a 2D dictionary that contains all the possible combinations of a read with a variant and the number of mismatches.
            for var in var_list:
                    for key in d:
                            temp = d[key]
                            if var in temp:
                                    index = temp.index(var)
                                    val = d_1[key][index]
                                    #only for perfect matches.
                                    '''if val >3:
                                    	x[key][var] = -1
                                    '''
                                    x[key][var] = int(val)

            df = pd.DataFrame(x).T.fillna(-1)
    return df

def tree():
    return defaultdict(tree)

#Compute P(read j | var i)
def compute_probability(n, k):
    b = comb(n, k, exact=False)
    x = math.pow(0.99,(n-k))
    y = math.pow(0.01,k)
    prob = x*y*b
    return prob

#Compute proportions of a variant
def compute_proportions(dataframe):
    prob_list = [] #a list to hold the probabilities
    for row in dataframe.itertuples(index=False):
            temp_list = list(row)
            #compute the probability for each row in the matrix
            for i in range(len(temp_list)):
                    if temp_list[i] >= 0:
                            temp_list[i] = compute_probability(76,int(temp_list[i]))
                    else:
                            temp_list[i] = 0
            total = sum(temp_list)
            #solve for k
            temp_list = [j*(100000000.0/(100000000*total)) for j in temp_list]
            prob_list.append(temp_list)
    col_sums = [sum(k) for k in zip(*prob_list)]
    total_sum = sum(col_sums)
    prop_list = [100.0*l*(1/total_sum) for l in col_sums]
    return prop_list    

def create_dictionary(keys, vals):
    my_dict = dict()
    if len(keys) == len(vals):
            for i in range(len(keys)):
                    my_dict[keys[i]] = vals[i]
    return my_dict 

def compute_likelihood(df):
    numVar = df.shape[1]
    likelihood_list = list()
    max_mm = 2
    
    for row in df.itertuples(index=False):
        read = list(row)
        
        temp = list()
        for i in range(numVar):
            if read[i] == -1:
                prob = (0.01)**(max_mm+1) * (0.99)**(76 - max_mm -1)
                temp.append(prob)
            else:
                prob = (0.01)**(read[i]) * (0.99)**(76 - read[i])
                temp.append(prob)
                
        likelihood_list.append( sum(temp) )
    
    likelihood_list = [i/(2.0*76*numVar) for i in likelihood_list]
    neg_log_likelihood = [-1.0*np.log10(j) for j in likelihood_list]
    
    score = sum(neg_log_likelihood)
    return score

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--path", required = True, help = "Path to data table")
ap.add_argument("-g", "--gene", required = True,  help="name of gene")
args = vars(ap.parse_args())

path = args["path"]

#generate matrix
df = Generate_Matrix(path)
df.rename(columns={'Unnamed: 0': 'Read'}, inplace=True)
#predict variants
pred_object_val,var_predicted,reads_cov, all_solutions, all_objective = ilp.solver(df)
#print var_predicted
#print all_solutions
df2 = df.loc[reads_cov,var_predicted]

#write matrix to file
df2.to_csv(args["gene"]+'_predicted_matrix.csv', sep='\t')
fig, ax = plt.subplots()
df2.hist(bins=np.histogram(df2.values.ravel())[1],ax=ax)
fig.savefig(args["gene"]+'_matrix_plots.png')

#compute proportions
prop = compute_proportions(df2)
pred_prop = create_dictionary(var_predicted, prop)

#score list and proportions
score_list = list()
min_score = sys.maxint
min_sol_list = list()
print("Total number of optimal solutions: {}".format(len(all_solutions)))

for i in range(len(all_solutions)):
    print("Solution: {}".format(all_solutions[i]))
    print("Objective value: {}".format(all_objective[i]))
    print("Proportion: {}".format(compute_proportions(df.loc[reads_cov, all_solutions[i]])))
    score = compute_likelihood(df.loc[reads_cov, all_solutions[i]])
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