#!/usr/bin/python
'''
This script take reads(text file) specific to a sample as input, gene name.
It outputs identified variants at that gene and their proportions.
'''
from __future__ import division
from collections import defaultdict
from scipy.misc import comb
import pandas as pd
import math
import variantILP as varSolver
import numpy as np
import getStrAndProp as gsp
import itertools
import cplex
import sys
import matplotlib.pyplot as plt
import re

NO_BINOM = False

def ConvertAscii(Qscore):
    result = []
    for char in Qscore:
        #subtract 33 for predefined offset
        result.append(ord(char))
    return result

def compute_QAvg(Qmatrix):
    Qavg = Qmatrix.sum()
    nonZero = Qmatrix.astype(bool).sum(axis=1).sum()
    Qavg = Qavg.sum()/nonZero
    return Qavg

def compute_QSum(Qmatrix):
    return (Qmatrix.sum()).sum()

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
    prob_list = [0.0]*dataframe.shape[1]
    
    for row in dataframe.itertuples(index=False):
        mmInfo = [i for i in list(row) if i>=0]
        min_mm = min(mmInfo)
        numOfVar_minMm = len([i for i in list(row) if i== min_mm])
        
        for i in range(len(list(row))):
            if list(row)[i] == min_mm:
                prob_list[i] += 1/numOfVar_minMm
                
    normalize_term = 1.0/(sum(prob_list))
    prob_list = [normalize_term * i for i in prob_list]
    return np.round(prob_list,10) 

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

def returnMismatchMatrix(path, option):
    df = pd.read_csv(path, sep='\t', header=None, usecols=[0,1,3], names=["Read", "Allele", "Mismatch"])
    df["Mismatch"] = df["Mismatch"].str[-1]     #grab mismatch number
    df["Mismatch"] = pd.to_numeric(df["Mismatch"], errors='coerce')
    
    if option == "paired":
        df = df.groupby(["Read", "Allele"], as_index=False)["Mismatch"].sum()   #combine two rows if option is paired
        
    matrix = df.pivot(index="Read", columns="Allele", values="Mismatch")    #transform into rows=reads, columns=alleles
    matrix = matrix.fillna(-1)
    
    if option == "paired":
        matrix = matrix[(matrix>= -1) & (matrix<=6)]
    else:
        matrix = matrix[(matrix>= -1) & (matrix<=3)]
    
    matrix = matrix.fillna(-1)
    matrix = matrix[(matrix.T != -1).any()]     #remove any rows with all -1 i.e. reads do not map to any alleles after limiting mm
    matrix = matrix.loc[:, (matrix != -1).any(axis=0)]  #remove any alleles not mapped by any reads after limiting mm
    return matrix

def writeReadTable(capGene, iteration, option):
    readOutFile = open("{0}_{1}_{2}NoHeader.sam".format(capGene, iteration, option))
    writefile = open("{0}_{1}_{2}_reads.txt".format(capGene, iteration, option), "w")
    for line in readOutFile:
        fields = line.strip("\t").split()
        read = fields[0]
        allele = fields[2]
        quality = fields[10]
        mm = [i for i in fields if i.startswith("XM:i:")][0]  #bowtie2
#        mm = [i for i in fields if i.startswith("NM")][0]   #bowtie
        mm_pos = [j for j in fields if j.startswith("MD:Z:")][0]
        
        writefile.write(read + "\t" + allele + "\t" + quality + "\t" + mm + "\t" + mm_pos + '\n')
        
    readOutFile.close()
    writefile.close()
    
def combiningTag(a, b):
    firstNumbers = re.split("\D+", a)
    firstNumbers = list(filter(None, firstNumbers))
    first = firstNumbers[-1]
    secondNumbers = re.split("\D+", b)
    secondNumbers = list(filter(None, secondNumbers))
    second = secondNumbers[0]
    numbChanged = int(first) + int(second)
    
    combined = (a)[:-len(firstNumbers[-1])] + str(numbChanged) + (b)[len(secondNumbers[0]):]
    return combined

def reconstructMDTag(md):
    if '^' in md:
        fields = re.split("\^[ATCG]+", md)
        fields = list(filter(None, fields))
        
        appendedStr = fields[0]
        for i in range(1, len(fields) ):
            appendedStr = combiningTag(appendedStr, fields[i])
            
    else:
        appendedStr = md
        
    return appendedStr

def returnQuality(quality, mm_pos):
    q_list = list()
    
    for index in range(len(mm_pos)):
        temp = re.split("\D+", mm_pos[index])
        temp = [int(i) for i in temp]
#        print(index)
#        print(mm_pos[index])
        
        calculate_quality_pos = list()
        calculate_quality_pos.append(temp[0])
        for j in range(1, len(temp)-1):
            calculate_quality_pos.append(calculate_quality_pos[j-1] + temp[j] + 1)
#        print(calculate_quality_pos)
        q = [ord( (quality[index])[k] ) for k in calculate_quality_pos]
        q_list.append(sum(q))
            
    return q_list

def returnQualityMatrix(path, option):
    df = pd.read_csv(path, sep='\t', header=None, usecols=[0,1,2,3,4], names=["Read", "Allele", "Quality", "Mismatch", "Mm position"])
    df.loc[:, "Mismatch"] = df.loc[:, "Mismatch"].str[-1]
    df.loc[:, "Mismatch"] = pd.to_numeric(df.loc[:, "Mismatch"], errors='coerce')
    df["Mm position"] = df["Mm position"].str.extract("MD:Z:(.*)", expand=False)
    zeroMismatch = (df["Mismatch"] == 0)
    df["Mm position"] = df["Mm position"].apply(reconstructMDTag)
    df.loc[~zeroMismatch, "Quality"] = returnQuality(df[~zeroMismatch]["Quality"].tolist(), df[~zeroMismatch]["Mm position"].tolist())
    df.loc[zeroMismatch, "Quality"] = 0
    df["Quality"] = pd.to_numeric(df["Quality"], errors='coerce')
    
    if option == "paired":
        tempDF = df.groupby(["Read", "Allele"], as_index=False)["Mismatch", "Quality"].sum()
        tempDF = tempDF[(tempDF["Mismatch"] >= 0) & (tempDF["Mismatch"]<=6)]
    else:
        tempDF = df[["Read", "Allele", "Mismatch", "Quality"]]
        tempDF = tempDF[(tempDF["Mismatch"] >= 0) & (tempDF["Mismatch"]<=3)]
    
    tempDF.reset_index(inplace=True, drop=True)    
    matrix = tempDF.pivot(index="Read", columns="Allele", values="Quality")
    
    if option == "paired":
        matrix = matrix.fillna(186)
    else:
        matrix = matrix.fillna(93)
        
    return matrix

#Only return solutions with minimum number of variants
def getVarAndProp(gene, paired_path, samp):
    #generate matrix
    pairedDF = returnMismatchMatrix(paired_path, "paired")
#    singletonDF = returnMismatchMatrix(singleton_path, "singleton")
#    dataMatrixDF = pd.concat([pairedDF, singletonDF])
#    dataMatrixDF = dataMatrixDF.fillna(-1)
    dataMatrixDF = pairedDF

    #Generate quality matrix
    Qmatrix_paired = returnQualityMatrix(paired_path, "paired")
    Qmatrix = Qmatrix_paired
#    Qmatrix_singleton = returnQualityMatrix(singleton_path, "singleton")
#    Qmatrix = pd.concat([Qmatrix_paired, Qmatrix_singleton])
#    Qmatrix.loc[Qmatrix_paired.index] = Qmatrix.loc[Qmatrix_paired.index].fillna(186)
#    Qmatrix.loc[Qmatrix_singleton.index] = Qmatrix.loc[Qmatrix_singleton.index].fillna(93)
    
    #predict variants
    pred_object_val,var_predicted,reads_cov, all_solutions, all_objective = varSolver.solver(dataMatrixDF)
    
    #Calculate quality scores of solutions
    Qscores = list()
#    minVar_Qscores = list()
    for i in range(len(all_solutions)):
        Qscores.append(compute_QSum(Qmatrix.loc[reads_cov,all_solutions[i]]))
        
    print("Quality scores for all solutions: {}".format(Qscores))
    print("Solutions: {}".format(all_solutions))
    
    #The required alleles
    dataMatrix_pred = dataMatrixDF.loc[reads_cov,var_predicted]
#    minVar_solutions = [sol for sol in all_solutions if len(sol) == min(map(len,all_solutions))]
    minVar_solutions = all_solutions
    
    #Calculate quality scores for minimum alleles
#    for i in range(len(minVar_solutions)):
#        minVar_Qscores.append(compute_QSum(Qmatrix.loc[reads_cov, minVar_solutions[i]]))
        
#    print("Quality scores for solutions with minimum alleles: {}".format(minVar_Qscores))
#    print("Minimum alleles solutions: {}".format(minVar_solutions))
    
    #Return solutions with minimum quality scores
    min_qscore = min(Qscores)
    minQscore_sol = [minVar_solutions[i] for i in range(len(minVar_solutions)) if Qscores[i] == min_qscore]
    print("Min quality score: {}".format(min_qscore))
    print("Solutions with minimum quality score: {}".format(minQscore_sol))
    if len(minQscore_sol) > 1:
        print("@@@@@@@@@@@@ More than 1 solution having minimum quality score @@@@@@@@@@@@")
    
    #Likelihood approach 
    ''' 
    #compute proportions
    prop = compute_proportions(dataMatrix_pred)
    pred_prop = create_dictionary(var_predicted, prop)
    
    #score list and proportions
    score_list = list()
    min_score = sys.maxint
    
    for i in range(len(minVar_solutions)):
        score = compute_likelihood(dataMatrixDF.loc[reads_cov, minVar_solutions[i]])
        score_list.append(score)
        
        if score <= min_score:
            min_score = score
            
    #Give some names to the solutions for further identifications and get the indices for sorted likelihood list
    #Sort quality score and produce some graphs
    sortedIndex_score_list = np.argsort(score_list)
    sortedIndex_qscore_minVarSol = np.argsort(Qscores)
    likelihood_score_dict = dict()
    qscores_minVarSol_dict = dict()
    sol_name_dict = dict()
    
    for i in range(len(minVar_solutions)):
        sol_name_dict["sol_{}".format(i)] = minVar_solutions[sortedIndex_score_list[i]]
        likelihood_score_dict["sol_{}".format(i)] = score_list[sortedIndex_score_list[i]]
        qscores_minVarSol_dict["sol_{}".format(i)] = Qscores[sortedIndex_qscore_minVarSol[i]]
            
#    plt.figure()
#    sorted_solution_namelist = ["sol_{}".format(i) for i in range(len(minVar_solutions))]
#    plt.xticks(range(len(minVar_solutions)), sorted_solution_namelist, rotation=20)
#    plt.scatter(range(len(minVar_solutions)), [likelihood_score_dict[name] for name in sorted_solution_namelist], s=50)
#    plt.xlabel('Solution i ')
#    plt.ylabel('Negative log likelihood')
#    plt.savefig("{0}_{1}_sol_likelihood".format(samp, gene))
    
#    plt.figure()
#    plt.xticks(range(len(minVar_solutions)), sorted_solution_namelist, rotation=20)
#    plt.scatter(range(len(minVar_solutions)), [qscores_minVarSol_dict[name] for name in sorted_solution_namelist], s=50)
#    plt.xlabel('Solution i ')
#    plt.ylabel('Average quality score')
#    plt.savefig("{0}_{1}_minVarSol_avgQscore".format(samp, gene))
    
#    minVar_minNegLog_solutions = [minVar_solutions[i] for i in range(len(minVar_solutions)) if min_score <= score_list[i] <= 1.01*min_score]
    '''
    
    ''' ====== '''
    #compute proportions
    #solutionsAndProp_dict is a dictionary in which the keys are just indices and values are dictionaries, with variant as key and proportion as value
    solutionsAndProp_dict = dict()
    track=0
    for sol in minQscore_sol:
        dataMatrix_pred = dataMatrixDF.loc[reads_cov, sol]
        prop = compute_proportions(dataMatrix_pred)
        pred_prop = create_dictionary(sol, prop)
        solutionsAndProp_dict[track] = pred_prop
        track += 1
        
#    plt.close('all')        
    return solutionsAndProp_dict
    
    
    
