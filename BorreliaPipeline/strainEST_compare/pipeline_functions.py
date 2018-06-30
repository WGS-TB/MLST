#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:33:50 2017

@author: stanleygan

Project: Illuminating the diversity of pathogenic bacteria Borrelia Burgdorferi in tick samples

"""
from __future__ import division
from collections import defaultdict, Counter
from scipy.misc import comb
from scipy.stats import entropy
from math import log
import pandas as pd
import csv
import os
import re
import itertools
import numpy as np
import cplex
import sys
import math
import variantILP as varSolver
#import matplotlib.pyplot as plt
#from scipy.spatial.distance import hamming

errorThres=0.1

''' Extend class as need to retrieve solution if gap does not converge after certain time'''
class TimeLimitCallback(cplex.callbacks.MIPInfoCallback):
    def __call__(self):
        if not self.aborted and self.has_incumbent():
            relGap = 100.0 * self.get_MIP_relative_gap()
            totalTimeUsed = self.get_time() - self.starttime
            if totalTimeUsed > self.timelimit and relGap < self.acceptablegap:
                print("Good enough solution at", totalTimeUsed, "sec., gap =",
                      relGap, "%, quitting.")
                self.aborted = True
                self.abort()

''' ====================================== Functions related to processing raw data ======================================================= '''
'''
Return the sum of all quality scores in Qmatrix
Input: Qmatrix, dataframe
'''
def compute_QSum(Qmatrix):
    return (Qmatrix.sum()).sum()

'''
Compute the probability of read of length n mapping to a variant with k mismatches using the binomial distribution
Input: n, integer
        k, integer
'''
def compute_probability(n, k):
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

def kallisto_proportions(alleles, kal_cmd, seq_dict, first_fa, second_fa):
    with open("temp_prop.fas", "w") as f:
        for a in alleles:
            f.write(">"+a+"\n")
            f.write(seq_dict[">"+a]+"\n")

    #kallisto index
    print os.getcwd()
    kal_idx_cmd = kal_cmd + ' index -i temp_prop.idx temp_prop.fas >/dev/null 2>&1'
    os.system(kal_idx_cmd)

    #Run kallisto quantifier
    kallisto_cmd = kal_cmd + ' quant -t 4 -i temp_prop.idx -o temp_prop {0} {1} >/dev/null 2>&1'.format(first_fa, second_fa)
    os.system(kallisto_cmd) 
    output_file = pd.read_csv(os.path.join(os.getcwd(), 'temp_prop', 'abundance.tsv'),sep='\t')
    DF = output_file.loc[:,['target_id','est_counts']]
    DF = DF[DF['est_counts'] != 0]
    DF['est_counts'] = (DF['est_counts']/float(DF['est_counts'].sum()))
    #DF = DF[DF['est_counts'] > 1.0]
    #DF['est_counts'] = (DF['est_counts']/DF['est_counts'].sum())
    var_predicted = DF['target_id'].tolist()
    props = DF['est_counts'].tolist()

    prop_dict = {var_predicted[i]:props[i] for i in range(len(var_predicted))}

    return prop_dict
    
    
'''
Input: Dataframe with rows=reads, columns=variants
Output: The proportions of variants (type list)
'''
def compute_proportions(dataframe):
    #computes the proportion of a set of variants given a set of reads uing probabilistic methods
    prob_list = [0.0]*dataframe.shape[1]
    for row in dataframe.itertuples(index=False):
        #mmInfo = [i for i in list(row) if i>=0]
        mmInfo = [i for i in list(row) if i!=6]
        min_mm = min(mmInfo)
        numOfVar_minMm = len([i for i in list(row) if i== min_mm])
        
        if numOfVar_minMm != len(list(row)) or len(list(row)) == 1:
        #if numOfVar_minMm == 1:
            for i in range(len(list(row))):
                if list(row)[i] == min_mm:
                    prob_list[i] += 1/numOfVar_minMm

    ##Only run entropy if there are more than 2 alleles
    #if mismatch_df.shape[1] != 1:            
    #    #compute entropy
    #    matrix_mismatch = mismatch_df.as_matrix()
    #    filtered_matrix_mismatch = list()
    #    for i in range(matrix_mismatch.shape[0]):
    #        temp = set(matrix_mismatch[i,:])
    #        if len(temp) > 1:
    #            filtered_matrix_mismatch.append(list(matrix_mismatch[i,:]))

    #    filtered_matrix_mismatch = np.array(filtered_matrix_mismatch)

    #    entropy_list = list()
    #    for i in range(filtered_matrix_mismatch.shape[1]):
    #        distrib = list(filtered_matrix_mismatch[:,i])        

    #        #In case we encounter -1, we insert different values for -1
    #        new_distrib = list()
    #        dummy = max(distrib) + 1
    #        for d in distrib:
    #            if d == -1:
    #                new_distrib.append(dummy)
    #                dummy += 1
    #            else:
    #                new_distrib.append(d)
    #    
    #        new_distrib_dict = dict(Counter(new_distrib))
    #        new_distrib_prop = [float(new_distrib_dict[k])/len(new_distrib) for k in new_distrib_dict.keys()]
    #        entro = entropy(new_distrib_prop, base=2)
    #        entropy_list.append(entro)
    #    
    #    scaled_entropy_list = [i/log(filtered_matrix_mismatch.shape[0], 2) for i in entropy_list]
    #    weights = [1-i for i in scaled_entropy_list]
    #else:
    #    weights=[1]*len(prob_list)

    ##multiply probs with weights based on entropy
    #prob_list = [w*p for (w, p) in itertools.izip(weights, prob_list)]
    
    try:            
        normalize_term = 1.0/(sum(prob_list))
    except ZeroDivisionError:
        print("no unique first score")
        normalized_term = 0
        
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
Input: A dataframe with rows=reads, columns=variants, max_mm=maximum mismatch set
Output: Negative log likelihood score of this solution
'''    
def compute_likelihood(df, max_mm):
    numVar = df.shape[1]
    likelihood_list = list()
    
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

'''
Return a mismatch dataframe where row=reads and columns=alleles, entries=# of mismatches
Input: path, absolute path to reads.txt file
        option, "paired" if combining mate pair reads into 1
'''
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

'''
Summarize required information from SAM file, where the reads.txt file contains information
for each read mapped, alleles that it map to, base quality score, number of mismatches and mismatch position
'''
#def writeReadTable(capGene, iteration, option):
#    readOutFile = open("{0}_{1}_{2}NoHeader.sam".format(capGene, iteration, option))
#    writefile = open("{0}_{1}_{2}_reads.txt".format(capGene, iteration, option), "w")
#    for line in readOutFile:
#        fields = line.strip("\t").split()
#        read = fields[0]
#        allele = fields[2]
#        quality = fields[10]
#        mm = [i for i in fields if i.startswith("XM:i:")][0]  #bowtie2
##        mm = [i for i in fields if i.startswith("NM")][0]   #bowtie
#        mm_pos = [j for j in fields if j.startswith("MD:Z:")][0]
#        
#        writefile.write(read + "\t" + allele + "\t" + quality + "\t" + mm + "\t" + mm_pos + '\n')
#        
#    readOutFile.close()
#    writefile.close()
   
'''
Combine two tags for mismatch position. For example, if mismatch position is "16^A2G2", we would like to get "18G2" as we do not care about
insertion details
Input: a, first string
        b, second string
'''
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

#Reconstruct the MD tag in SAM file as we are only interested in mismatches information but not insertion
#Input: md, md tag which is a string
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

#Return the base quality according to the position specified in MD tag
#Input: quality, a string
#       mm_pos, a string
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
        q = [(k-33)/93 for k in q]
        q_list.append(sum(q))

    return q_list

'''
Return a dataframe where rows=reads and columns=alleles, entries=quality score
Input: path, absolute path to the reads.txt file
        option, "paired" if combining mate pairs read as 1
'''
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
    #The max quality score is 93.As we limit to 3 mismatches, hence the maximum of an entry is 93*3
    if option == "paired":
        matrix = matrix.fillna(6)
    else:
        matrix = matrix.fillna(3)
        
    return matrix

''' ====================================== Functions related to allele prediction ======================================================= '''

'''
Return a dictionary with key=indices, values=dictionary where key=allele, value=proportions
Input:
    gene, name of gene
    paired_path, absolute path to reads.txt file
    samp, sample name
'''
def getVarAndProp(gene, paired_path, samp):
    #generate matrix
    dataMatrixDF = returnMismatchMatrix(paired_path, "paired")

    #Generate quality matrix
    Qmatrix = returnQualityMatrix(paired_path, "paired")
    
    #predict variants
    pred_object_val,var_predicted,reads_cov, all_solutions, all_objective = varSolver.solver(dataMatrixDF, Qmatrix, "paired")
    
    score_list = list()
    min_score = sys.maxint

    #Compute negative log likelihood score for each solution
    for i in range(len(all_solutions)):
        score = compute_likelihood(dataMatrixDF.loc[reads_cov, all_solutions[i]], 6)
        score_list.append(score)


        if score <= min_score:
            min_score = score

    argmin_score_list = [i for i in range(len(all_solutions)) if score_list[i] == min_score]
    if len(argmin_score_list) > 1:
        print("More than 1 solution having minimum negative log likelihood score.")
        lexico_min_score_sol = [all_solutions[i] for i in argmin_score_list]
        lexico_min_score_sol = sorted(lexico_min_score_sol)
        var_predicted = lexico_min_score_sol[0]
    else:
        var_predicted = all_solutions[argmin_score_list[0]]

    
    ''' ====== '''
    #compute proportions
    #solutionsAndProp_dict is a dictionary in which the keys are just indices and values are dictionaries, with variant as key and proportion as value
    solutionsAndProp_dict = dict()
    dataMatrix_pred = Qmatrix.loc[reads_cov, var_predicted]
    prop = compute_proportions(dataMatrix_pred)
    pred_prop = create_dictionary(var_predicted, prop)
    solutionsAndProp_dict[0] = pred_prop
    
    print("Solutions:{}".format(all_solutions))
    print("Score:{}".format(score_list)) 
    return solutionsAndProp_dict, len(all_solutions)

'''
This function is needed only if compatibility and quality score filtering are not discriminative enough. 
localMILP returns the objective value of the MILP for a given distribution
Input:
    sample, name of sample
    loci, a list of loci
    gene_solProp_dict, a dictionary with key=gene, values=dictionary where key=allele, value=proportion
    reference, a dataframe of all the existing strains
    objectiveOption: "all" means include all objective components, "noPropAndErr" omits proportion and error terms
'''
def localMILP(sample, loci, gene_solProp_dict, reference, objectiveOption, timelimit, gap):
    genesDF = pd.DataFrame(columns=loci)
    
    for gene in loci:
        genesDF[gene] = [gene_solProp_dict[gene]]
    
    data = dict()
    data[sample] = genesDF
    data = roundProp(data)
    newNameToOriName = dict()
    namingIndex=1
    for i in sorted(data.keys()):
        newNameToOriName["s{}".format(namingIndex)] = i
        data["s{}".format(namingIndex)] = data.pop(i) 
        namingIndex += 1
    
    ''' ============================================== Data handling ====================================================== '''
    #paramaters
    propFormat = 1    #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    allSamples = data.keys()
        
    #Get proportions of variants at different locus for each sample
    varAndProp = returnVarAndProportions(data)
    
    #Get the combinations at all loci across all samples
    strains, numOfComb= returnCombinationsAndNumComb(data, numLoci, loci)
#    numOfComb = strainAndNumComb[1]
    uniqueStrains = strains.drop_duplicates(loci)
    uniqueStrains = (uniqueStrains[loci]).reset_index(drop=True)
    uniqueStrains["ST"] = uniqueStrains.index.values + 1    #assign indices for strains or each unique combinations
    strains = strains.merge(uniqueStrains, indicator=True, how="left")    #assign the index to the combinations(as strain data frame contains duplicated rows)
    strains = strains.drop("_merge",1)
    
    #For each variants, get a mapping of which strains it maps to
    varSampToST = mapVarAndSampleToStrain(strains, loci, allSamples)
    
    #weights and decision variables for proportion of strains. weight=0 if the strain is in reference, otherwise =1. Notice there will be duplications of strain types
    #here because for proportions, we consider sample by sample rather than unique strain types
    proportionWeightDecVarDF = strains.merge(reference, indicator=True, how="left")
    proportionWeightDecVarDF["_merge"].replace(to_replace="both", value=0, inplace=True)
    proportionWeightDecVarDF["_merge"].replace(to_replace="left_only", value=1, inplace=True)
    proportionWeightDecVarDF = proportionWeightDecVarDF.rename(columns = {"_merge":"Weights"})
    
    #Add proportion decision variable names
    proportionWeightDecVarDF["Decision Variable"] = np.nan
    
    for samp in allSamples:
        thisSample = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Sample']
        propNameTemp = ["pi_%s_%d" %t for t in itertools.izip(thisSample, range(1,1+thisSample.shape[0]))]
        #shorter name as CPLEX can't hold name with >16 char. Use last 3 digits of sample name to name decision variables i.e. SRR2034333 -> use 333
        propNameTemp = [ele.replace("pi_{}".format(samp), "pi_{}".format(samp)) for ele in propNameTemp]   
        proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp, 'Decision Variable'] = propNameTemp
        
    #weights and decision variables for unique strain types, weight=0 if strain is in reference, otherwise=1. no duplications
    strainWeightDecVarDF = proportionWeightDecVarDF.drop_duplicates(loci)
    retainCol = loci + ['Weights', 'ST']
    strainWeightDecVarDF = strainWeightDecVarDF[retainCol].reset_index(drop=True)
    strainWeightDecVarDF["Decision Variable"] = ["a{}".format(i) for i in range(1, strainWeightDecVarDF.shape[0] + 1)]
    
    '''==================================== Forming ILP here ================================================'''
     #Form a CPLEX model
    model = cplex.Cplex()
    
    #Some bound on cplex solver when gap finds it hard to converge
    timelim_cb = model.register_callback(TimeLimitCallback)
    timelim_cb.starttime = model.get_time()
    timelim_cb.timelimit = timelimit
    timelim_cb.acceptablegap = gap
    timelim_cb.aborted = False
    
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=strainWeightDecVarDF['Weights'].values.tolist(), names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
    #add proportions decision variables
    if objectiveOption == "noPropAndErr":
        model.variables.add(lb=[0]*proportionWeightDecVarDF.shape[0], ub=[propFormat]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))
    else:
        model.variables.add(obj=[i for i in proportionWeightDecVarDF['Weights'].values.tolist()],lb=[0]*proportionWeightDecVarDF.shape[0],ub=[propFormat]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))
        
    #add linear constraints such that for each sample, the sum of the proportions of its variants combination = 1
    propVarSumTo1 = list()
    
    for samp in allSamples:
        temp = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Decision Variable'].tolist()        
        propVarSumTo1.append([temp, [1]* len(temp)])
        
    model.linear_constraints.add(lin_expr=propVarSumTo1, rhs=[propFormat]*len(propVarSumTo1), senses=["E"]*len(propVarSumTo1), names=["c{0}".format(i+1) for i in range(len(propVarSumTo1))])                                                                 
    
    #add linear constraints such that each decision variable a_i must be at least pi_jk in which pi_jk is the proportion of V_jk and V_jk=a_i
    #By this, if we use any of the pi, we force a_i to be 1
    indicLargerPropDF = pd.DataFrame(columns=["ST","Indicator"])
    indicLargerPropDF["ST"] = strainWeightDecVarDF["ST"]
    indicLargerPropDF["Indicator"] = strainWeightDecVarDF["Decision Variable"]
    indicLargerPropDF = (indicLargerPropDF.merge(proportionWeightDecVarDF, indicator=True, how="left", on="ST"))[["ST","Indicator","Decision Variable"]]
    indicLargerPropDF.rename(columns={"Decision Variable": "Proportion Variable"}, inplace=True)
    indicMinusProp = list()
    for i,pi in itertools.izip(indicLargerPropDF["Indicator"].tolist(), indicLargerPropDF["Proportion Variable"].tolist()):
        indicMinusProp.append([[i, pi],[propFormat, -1]])  
    
    model.linear_constraints.add(lin_expr=indicMinusProp, rhs=[0]*len(indicMinusProp), senses=["G"]*len(indicMinusProp), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusProp))] )
    
    #Also, add linear constraints such that a_i - average of pi_jk <= 0.999. Otherwise will have case that a_i=1 and for all pi_jk, pi_jk=0
    indicMinusAvgPropLess1_DF = indicLargerPropDF.groupby("Indicator")["Proportion Variable"].apply(list).reset_index()
    indic = indicMinusAvgPropLess1_DF["Indicator"].tolist()
    pV = indicMinusAvgPropLess1_DF["Proportion Variable"].tolist()
    indicMinusAvgPropLess1_LHS = list()
    
    for i in range(len(indic)):
        a_i = indic[i]
        pi_i = pV[i]
        temp = list()
        size = len(pi_i)
        temp.append(a_i)
        coef = list()
        coef.append(propFormat)
        
        for j in range(size):
            temp.append(pi_i[j])
            coef.append(-1.0/size)
            
        indicMinusAvgPropLess1_LHS.append([temp, coef])
    
    tolerance = 0.01*propFormat*0.01     #how much tolerance we set for the upper bound    
    model.linear_constraints.add(lin_expr=indicMinusAvgPropLess1_LHS, rhs=[propFormat - tolerance]*len(indicMinusAvgPropLess1_LHS), senses=["L"]*len(indicMinusAvgPropLess1_LHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusAvgPropLess1_LHS))])
    model.linear_constraints.add(lin_expr=indicMinusAvgPropLess1_LHS, rhs=[0]*len(indicMinusAvgPropLess1_LHS), senses=["G"]*len(indicMinusAvgPropLess1_LHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusAvgPropLess1_LHS))])
    
    #add error variables and linear constraints related to error terms
    #create error variable names
    varAndProp["Decision Variable"] = ["d_{}_".format(samp) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Decision Variable"] = varAndProp["Decision Variable"] + varAndProp["Variant"]
    
    #add error variables
    #errorThres = 0.1
    model.variables.add(obj=[1]*varAndProp.shape[0], names=varAndProp["Decision Variable"].tolist(), lb= [0]*varAndProp.shape[0], ub= [errorThres]*varAndProp.shape[0], types=[model.variables.type.continuous]*varAndProp.shape[0])

    #add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
    piDotComb = list()
    piDotComb_2 = list()
    propConstrRHS = list()
    for locusName in varSampToST:
        temp=list()
        varSampToSTDict = varSampToST[locusName][0]
        for (var, sample) in varSampToSTDict:
            strainTypes = varSampToSTDict[(var, sample)]
            propDecVar = proportionWeightDecVarDF[(proportionWeightDecVarDF["ST"].isin(strainTypes)) & (proportionWeightDecVarDF["Sample"] == "{}".format(sample))]["Decision Variable"]
            errorDecVar = varAndProp[(varAndProp["Variant"] == var) & (varAndProp["Sample"] == sample)]["Decision Variable"]
            propConstrRHS.append(  float( ( (data["{}".format(sample)])[locusName][0] )[var] )  )
            piDotComb.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [-1]])
            piDotComb_2.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [1]])
           
    model.linear_constraints.add(lin_expr=piDotComb, rhs=propConstrRHS, senses=["L"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])     
    model.linear_constraints.add(lin_expr=piDotComb_2, rhs=propConstrRHS, senses=["G"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))]) 

    #Export some info for MATLAB use
    #writeInfoToCsv()
    
    ''' ================================== Solve ILP ========================================== '''
    #model.write("borreliaLP.lp")
#    model.set_results_stream(None)
    model.solve()
#    model.write("{}.lp".format(sample))
    
    #options for searching more optimal solutions
#    model.parameters.mip.pool.capacity.set(50)
#    model.parameters.mip.pool.intensity.set(4)
#    model.parameters.mip.limits.populate.set(100)
#    model.parameters.mip.pool.absgap.set(0)
#    model.parameters.mip.pool.replace.set(1)
#    model.populate_solution_pool()
    
    objvalue = model.solution.get_objective_value()
    varNames = model.variables.get_names()
    varValues = model.solution.get_values(varNames)
    conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
    conclusion["Decision Variable"] = varNames
    conclusion["Value"] = varValues
    error = conclusion[conclusion["Decision Variable"].str.contains("d_")]
    nonZero_error = error[(error["Value"] <= -0.001) | (error["Value"] >= 0.001)]
    if nonZero_error.shape[0] != 0:
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^ error ^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print nonZero_error
        
#    objStr = conclusion[conclusion["Decision Variable"].str.contains("^a")]["Decision Variable"].tolist()
#    objProp = conclusion[conclusion["Decision Variable"].str.contains("^pi")]["Decision Variable"].tolist()
#    objErr = conclusion[conclusion["Decision Variable"].str.contains("^d")]["Decision Variable"].tolist()
#    
#    objStr_coeff = model.objective.get_linear(objStr)
#    objProp_coeff = model.objective.get_linear(objProp)
#    objErr_coeff = model.objective.get_linear(objErr)
#    
#    sum_str = sum([val*coeff for val, coeff in itertools.izip(model.solution.get_values(objStr), objStr_coeff)])
#    sum_prop = sum([val*coeff for val, coeff in itertools.izip(model.solution.get_values(objProp), objProp_coeff)] )
#    sum_err = sum([val*coeff for val, coeff in itertools.izip(model.solution.get_values(objErr), objErr_coeff)] )
#    print("Objective value: {}".format(objvalue))
#    print("Strain component: {}".format(sum_str))
#    print("Prop component: {}".format(sum_prop))
#    print("Error component: {}".format(sum_err))
#    print("Sum :{}".format(sum_str+sum_prop+sum_err))
    
    return objvalue

'''
This function is needed only if compatibility and quality score filtering are not discriminative enough. 
localILP returns the objective value of the ILP for a given distribution
This ILP only consider strains, not taking into consideration any proportions involved
Input:
    sample, name of sample
    loci, a list of loci
    gene_solProp_dict, a dictionary with key=gene, values=dictionary where key=allele, value=proportion
    reference, a dataframe of all the existing strains
Output:
    solution_dict, dictionary where key=indices, value=dataframe related to that solution(information such as alleles at each locus)
    objective_dict, dictionary where key=indices, value=objective value of the i-th solution
    data, dictionary where key=sample name, value=dataframe which contains information about alleles and proportion at each locus
    strains, dataframe of unique strains for later use of localLP
'''
def localILP(sample, loci, gene_solProp_dict, reference):
    genesDF = pd.DataFrame(columns=loci)
    
    for gene in loci:
        genesDF[gene] = [gene_solProp_dict[gene]]
    
    data = dict()
    data[sample] = genesDF
    data = roundProp(data)
    newNameToOriName = dict()
    namingIndex=1
    for i in sorted(data.keys()):
        newNameToOriName["s{}".format(namingIndex)] = i
        data["s{}".format(namingIndex)] = data.pop(i) 
        namingIndex += 1
    allSamples = data.keys()
     #paramaters
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #Get the combinations at all loci across all samples
    strains, numOfComb = returnCombinationsAndNumComb(data, numLoci, loci)
    uniqueStrains = strains.drop_duplicates(loci)
    uniqueStrains = (uniqueStrains[loci]).reset_index(drop=True)
    uniqueStrains["ST"] = uniqueStrains.index.values + 1    #assign indices for strains or each unique combinations
    strains = strains.merge(uniqueStrains, indicator=True, how="left")    #assign the index to the combinations(as strain data frame contains duplicated rows)
    strains = strains.drop("_merge",1)
    
    #weights and decision variables for proportion of strains. weight=0 if the strain is in reference, otherwise =1. Notice there will be duplications of strain types
    #here because for proportions, we consider sample by sample rather than unique strain types
    strainWeightDecVarDF = strains.merge(reference, indicator=True, how="left")
    strainWeightDecVarDF["_merge"].replace(to_replace="both", value=0, inplace=True)
    strainWeightDecVarDF["_merge"].replace(to_replace="left_only", value=1, inplace=True)
    strainWeightDecVarDF = strainWeightDecVarDF.rename(columns = {"_merge":"Weights"})
    strainWeightDecVarDF = strainWeightDecVarDF.drop_duplicates(loci)
    retainCol = loci + ['Weights', 'ST']
    strainWeightDecVarDF = strainWeightDecVarDF[retainCol].reset_index(drop=True)
    strainWeightDecVarDF["Decision Variable"] = ["a{}".format(i) for i in range(1, strainWeightDecVarDF.shape[0] + 1)]
    
    #Relate sample and strain decision variable 
    samp_decVar_DF = strains.merge(strainWeightDecVarDF, how="left")[loci+["Sample", "Decision Variable", "ST"]]
    
    #For each allele, get a mapping of which strains it maps to
    varSampToST = mapVarAndSampleToStrain(samp_decVar_DF, loci, allSamples)
    
    '''==================================== Forming ILP here ================================================'''
    #Form a CPLEX model
    model = cplex.Cplex()
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=strainWeightDecVarDF['Weights'].values.tolist(), names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
    
    #Add linear constraints where strains chosen are able to describe all alleles seen in all samples
    #Add linear constraints where strains chosen are able to describe all alleles seen in all samples
    descAllAlleleLHS = list()
    for locusName in varSampToST:
        varSampToSTDict = varSampToST[locusName][0]
        for (var, sample) in varSampToSTDict:
            strainTypes = varSampToSTDict[(var, sample)]
            strainDecVar = samp_decVar_DF[(samp_decVar_DF["ST"].isin(strainTypes)) & (samp_decVar_DF["Sample"] == "{}".format(sample))]["Decision Variable"].tolist()
            descAllAlleleLHS.append([strainDecVar, [1]*len(strainDecVar)])
    
    model.linear_constraints.add(lin_expr=descAllAlleleLHS, rhs=[1]*len(descAllAlleleLHS), senses=["G"]*len(descAllAlleleLHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(descAllAlleleLHS))])

#    model.solve()
    
    #options for searching more optimal solutions
    #model.parameters.mip.pool.capacity.set(10)
#    model.set_results_stream(None)
    model.parameters.mip.pool.intensity.set(4)
#    model.parameters.mip.limits.populate.set(50)
    model.parameters.mip.pool.absgap.set(0)
    model.parameters.mip.pool.replace.set(1)
    model.populate_solution_pool()
    
    solution_dict = dict()
    objective_dict = dict()
    for i in range(model.solution.pool.get_num()):
        objvalue = model.solution.pool.get_objective_value(i)
        objective_dict[i] = objvalue
        varNames = model.variables.get_names()
        varValues = model.solution.pool.get_values(i,varNames)
        conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
        conclusion["Decision Variable"] = varNames
        conclusion["Value"] = varValues
    
        strainInfo = conclusion.merge(strainWeightDecVarDF[strainWeightDecVarDF["Decision Variable"].isin(varNames)])
        strainInfo["New/Existing"] = ["Existing" if w==0 else "New" for w in strainInfo["Weights"].tolist()]
        strainsNeeded = (strainInfo[strainInfo["Value"] == 1][loci + ["ST", "Weights"]])
        strainsNeeded.reset_index(drop=True, inplace=True)
        solution_dict[i] = strainsNeeded
        
#    print("Objective value: {}".format(objective_value))
    return solution_dict, objective_dict, data, strains, newNameToOriName

'''
This function is needed only if compatibility and quality score filtering are not discriminative enough. 
localLP returns the objective value of the LP for a given solution and the distribution.
This function takes solution from localILP and consider the effect of proportions
Input:
    solution, dataframe which contains alleles at each locus for a solution
    data, see localILP
    strains, see localILP
    reference, dataframe of existing strains
    loci, a list of locus
Output:
    objvalue, objective value for this solution in this LP
    feasible, indicator whether this solution is feasible
'''
def localLP(solution, data, strains, reference, loci, newNameToOriName):
    #paramaters
    propFormat = 1    #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    lociNames = list(reference.columns.values)
    numReference = reference.shape[0]
    allSamples = data.keys()
        
    #Get proportions of variants at different locus for each sample
    varAndProp = returnVarAndProportions(data)
    
    #Add propportion variables
    proportionWeightDecVarDF = strains.merge(solution, how='left', indicator=True)
    proportionWeightDecVarDF = proportionWeightDecVarDF[proportionWeightDecVarDF["_merge"] == "both"]
    proportionWeightDecVarDF.drop(["_merge"], axis=1, inplace=True)
    proportionWeightDecVarDF.reset_index(drop=True, inplace=True)
    
    #For each variants, get a mapping of which strains it maps to. Only consider those strains in given solution
    varSampToST = mapVarAndSampleToStrain(proportionWeightDecVarDF[loci+["Sample", "ST"]], loci, allSamples)
    
    #Add proportion variables names
    for samp in allSamples:
        thisSample = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Sample']
        propNameTemp = ["pi_%s_%d" %t for t in itertools.izip(thisSample, range(1,1+thisSample.shape[0]))]
        #shorter name as CPLEX can't hold name with >16 char. Use last 3 digits of sample name to name decision variables i.e. SRR2034333 -> use 333
        propNameTemp = [ele.replace("pi_{}".format(samp), "pi_{}".format(samp)) for ele in propNameTemp]   
        proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp, 'Decision Variable'] = propNameTemp
        
    ''' ===================================== Forming LP here =================================================== '''
 #Form a CPLEX model
    model = cplex.Cplex()
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=proportionWeightDecVarDF['Weights'].values.tolist(), lb=[0]*proportionWeightDecVarDF.shape[0], ub=[propFormat]*proportionWeightDecVarDF.shape[0], names=proportionWeightDecVarDF['Decision Variable'], types = [model.variables.type.continuous]* len(proportionWeightDecVarDF['Weights'].values.tolist()))
        
    #add linear constraints such that for each sample, the sum of the proportions of its variants combination = 1
    propVarSumTo1 = list()
    
    for samp in allSamples:
        temp = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Decision Variable'].tolist()        
        propVarSumTo1.append([temp, [1]* len(temp)])
        
    model.linear_constraints.add(lin_expr=propVarSumTo1, rhs=[propFormat]*len(propVarSumTo1), senses=["E"]*len(propVarSumTo1), names=["c{0}".format(i+1) for i in range(len(propVarSumTo1))])
    
    #add error variables and linear constraints related to error terms
    #create error variable names
    varAndProp["Decision Variable"] = ["d_{}_".format(samp) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Decision Variable"] = varAndProp["Decision Variable"] + varAndProp["Variant"]
    
    #add error variables
    #errorThres = 0.10
    model.variables.add(obj=[1]*varAndProp.shape[0], names=varAndProp["Decision Variable"].tolist(), lb=[0]*varAndProp.shape[0], ub= [errorThres]*varAndProp.shape[0], types=[model.variables.type.continuous]*varAndProp.shape[0])
    
   #add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
    piDotComb = list()
    piDotComb_2 = list()
    propConstrRHS = list()
    for locusName in varSampToST:
        temp=list()
        varSampToSTDict = varSampToST[locusName][0]
        for (var, sample) in varSampToSTDict:
            strainTypes = varSampToSTDict[(var, sample)]
            propDecVar = proportionWeightDecVarDF[(proportionWeightDecVarDF["ST"].isin(strainTypes)) & (proportionWeightDecVarDF["Sample"] == "{}".format(sample))]["Decision Variable"]
            errorDecVar = varAndProp[(varAndProp["Variant"] == var) & (varAndProp["Sample"] == sample)]["Decision Variable"]
            propConstrRHS.append(  float( ( (data["{}".format(sample)])[locusName][0] )[var] )  )
            piDotComb.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [-1]])
            piDotComb_2.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [1]])
           
    model.linear_constraints.add(lin_expr=piDotComb, rhs=propConstrRHS, senses=["L"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])           
    model.linear_constraints.add(lin_expr=piDotComb_2, rhs=propConstrRHS, senses=["G"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])                                                                   
    
    ''' ==== Solve ==== '''
    model.set_problem_type(0)   #set to LP problem
#    model.set_results_stream(None)
#    model.set_error_stream(None)
#    model.write("a.lp")
    model.solve()
#    print model.solution.get_status_string()
    feasible = False
    if model.solution.get_status() == 1:
        objvalue = model.solution.get_objective_value()
        feasible = True
    else:
        objvalue= -1
    
    objvalue = model.solution.get_objective_value()

    return objvalue, feasible

'''
Return the distribution which optimizes the local MILP (If more than 1, choose the one which is returned first)
Input:
    samp, sample name
    aTuple, tuple representing which distribution to consider
    aDict, a dictionary where key=gene, value=a dictionary where key=solution indices, value=a dictionary where key=allele, value=proportion of the allele
    This is confusing but here is an example: aDict = {clpA: {0:{clpA_1:0.5, clpA_2:0.5}, 1:{clpA_1:1.0}, 1:{...} }
                                                       clpX: {...}, ...}
    loci, a list of locus
    reference, dataframe of existing strains
    option, "all" if all objective components, "noPropAndErr" if omit proportion and error terms
Output:
    comb_minObjVal_dict, a dictionary where key=gene, value=a dictionary where key=allele, value=proportion
'''
def localMinimizer(samp, aTuple, aDict, loci, reference, option, timelimit, gap):
    track = 1
    objValue_list = list()
    print("\nNumber of combinations to run: {}\n".format(len(aTuple)))
    for combin in aTuple:
        print("\nxxxxxxxxxxxxxxxxx Combination : {} xxxxxxxxxxxxxxxxxxxxxxxxxxxx\n".format(track))
        comb_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, combin)}
        objVal = localMILP(samp, loci, comb_dict, reference, option, timelimit, gap)
        objValue_list.append(objVal)
        track += 1
        
    print("Objective Value: {}".format(objValue_list))
    #Choose the combination which has the lowest objective value
    minObjValIndex_list = np.argwhere(objValue_list == np.amin(objValue_list))
    minObjValIndex_list = minObjValIndex_list.flatten().tolist()
    if len(minObjValIndex_list) > 1:
        print("@@@@@@@@@@@@@@@@@@@@@@@ You have more than 1 distribution having same objective value @@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    
    minObjValIndex = minObjValIndex_list[0]
    comb_minObjVal = aTuple[minObjValIndex]
    comb_minObjVal_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, comb_minObjVal)}
    
    return comb_minObjVal_dict

'''
Heuristic: Return the distribution which optimizes the local ILP first, then the local LP (If more than 1, choose the one which is returned first)
Input:
    samp, sample name
    aTuple, tuple representing which distribution to consider
    aDict, a dictionary where key=gene, value=a dictionary where key=solution indices, value=a dictionary where key=allele, value=proportion of the allele
    This is confusing but here is an example: aDict = {clpA: {0:{clpA_1:0.5, clpA_2:0.5}, 1:{clpA_1:1.0}, 1:{...} }
                                                       clpX: {...}, ...}
    loci, a list of locus
    reference, dataframe of existing strains
Output:
    comb_minObjVal_dict, a dictionary where key=gene, value=a dictionary where key=allele, value=proportion
'''
def localMinimizer_sep(samp, aTuple, aDict, loci, reference):
    track = 1
    objValue_list = list()
    checkAllComb_feasibility = False
    print("\nNumber of combinations to run: {}\n".format(len(aTuple)))
    for combin in aTuple:
        print("\nxxxxxxxxxxxxxxxxx Combination : {} xxxxxxxxxxxxxxxxxxxxxxxxxxxx\n".format(track))
        comb_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, combin)}
        solution_dict, ilp_objective_dict, data, strains,newNameToOriName = localILP(samp, loci, comb_dict, reference)
        feasible_sol = list()
        lp_objective_dict = dict()
#        print solution_dict
        infeasibility = 0
        for i in solution_dict.keys():
            try:
                objvalue, feasible = localLP(solution_dict[i], data, strains, reference, loci, newNameToOriName)
                if feasible == False:
                    infeasibility += 1
                else:
                    feasible_sol.append(i)
                    lp_objective_dict[i] = objvalue
            except cplex.exceptions.errors.CplexSolverError as e:
                infeasibility += 1
        
        if infeasibility == len(solution_dict):
            print ("This combination has no feasible solutions")
        else:
            min_obj = np.inf
            for j in feasible_sol:
                if (ilp_objective_dict[j] + lp_objective_dict[j]) < min_obj:
                    min_obj = ilp_objective_dict[j] + lp_objective_dict[j]
            
            objValue_list.append(min_obj)
            print("Objective value: {}".format(min_obj))
            track += 1
            checkAllComb_feasibility = True
    
    if checkAllComb_feasibility == False:
        return -1
    
    print("Objective Value: {}".format(objValue_list))
    #Choose the combination which has the lowest objective value
    minObjValIndex_list = np.argwhere(objValue_list == np.amin(objValue_list))
    minObjValIndex_list = minObjValIndex_list.flatten().tolist()
    if len(minObjValIndex_list) > 1:
        print("@@@@@@@@@@@@@@@@@@@@@@@ You have more than 1 distribution having same objective value @@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    
    minObjValIndex = minObjValIndex_list[0]
    comb_minObjVal = aTuple[minObjValIndex]
    comb_minObjVal_dict = {gene: aDict[gene][i] for (gene, i) in itertools.izip(loci, comb_minObjVal)}
    
    return comb_minObjVal_dict

'''
Return tuples which are explainable by the existing strains
Input: 
    aTuple, a tuple of integer which represents which solution to consider at each locus
    gene_solProp_dict, see localMinimizer_sep
    loci, a list of locus
    reference, dataframe of existing strains
Output:
    compatible_tuples, a list of compatible tuples
'''
def compatibleFilter(aTuple, gene_solProp_dict, loci, reference):
    compatible_tuples = list()
    for combin in aTuple:
        comb_dict = {gene: gene_solProp_dict[gene][i] for (gene, i) in itertools.izip(loci, combin)}
        
        temp_boolean = True
        for allele in loci:
            temp_boolean = temp_boolean & reference[allele].isin(comb_dict[allele].keys())
            if sum(temp_boolean) == 0:
                break
            
        if sum(temp_boolean) != 0:
            compatible_tuples.append(combin)
            
    return compatible_tuples

''' ====================================== Functions related to strain prediction ======================================================= '''

'''Create a function to read all data files and return a dictionary where keys are the sample names, 
values are dataframes which contain the information about variants and proportions at different loci.
It also return total number of samples and starting sample number(based on last 3 digits)

Input
dataFilePath: File path that contains all your sample folders. dataFilePath should only contain directories of samples and 
              a reference csv
lociOrder: a list that contains order of columns that you want 
option: "all" if running on all samples
'''
def readData(dataFilePath, lociOrder, option):
    data = dict()
    sampleFold = list()
    
    if option == "all":     #dataFilePath = ...../variantsAndProp
        for folder in os.listdir(dataFilePath):
            if not folder.endswith(".csv"):
                sampleFold.append(folder)
    else:   #dataFilePath = ..../variantsAndProp/SRR2034333
        sampleFold = [dataFilePath.split("/")[-1]]
            
    numSamples = len(sampleFold)
#    startingSamp = min([int(i[-3:]) for i in sampleFold])
    
    for folder in sampleFold:
        data["{}".format(folder)] = pd.DataFrame(columns=lociOrder)   #require column to be a specfic order based on lociOrder
        if option == "all":
            sampleFilePath = dataFilePath + "/{}".format(folder)
        else:
            sampleFilePath = dataFilePath
        dirs= os.listdir(sampleFilePath)
        csvFiles = [i for i in dirs if i.endswith("proportions.csv")]
        
        temp = re.compile('(.*)_proportions.csv')  #grab gene name
        for f in csvFiles:
            gene = temp.findall(f)
            reader = csv.reader(open(sampleFilePath+"/"+f, 'r'))
            
            #store variants and respective proportions in a dictionary
            d = dict()  
            for variant, proportion in reader:
                d[variant] = proportion
            
            #wrap dictionary with a list for storage in dataframe
            alist =[d]
            (data["{}".format(folder)])[gene[0]] = alist
            
    return data, numSamples

#def readDataWithoutProp(dataFilePath, lociOrder, option):
#    data = dict()
#    sampleFold = list()
#    
#    if option == "all":     #dataFilePath = ...../variantsAndProp
#        for folder in os.listdir(dataFilePath):
#            if not folder.endswith(".csv"):
#                sampleFold.append(folder)
#    else:   #dataFilePath = ..../variantsAndProp/SRR2034333
#        sampleFold = [dataFilePath.split("/")[-1]]
#    
#    numSamples = len(sampleFold)
#    startingSamp = min([int(i[-3:]) for i in sampleFold])
#    
#    #data is a dictionary in which key=sample, value=a dataframe which has entries=list and the columns are loci
#    for folder in sampleFold:
#        data[folder] = pd.DataFrame(columns=lociOrder)
#        if option == "all":
#            sampleFilePath = dataFilePath + "/{}".format(folder)
#        else:
#            sampleFilePath = dataFilePath
#            
#        dirs= os.listdir(sampleFilePath)
#        csvFiles = [i for i in dirs if i.endswith("proportions.csv")]
#        temp = re.compile('(.*)_proportions.csv')  #grab gene name
#        for f in csvFiles:
#            gene = temp.findall(f)
#            reader = csv.reader(open(sampleFilePath+"/"+f, 'r'))
#            
#            #store alleles in a list
#            allele_list = list()
#            for row in reader:
#                allele_list.append(row[0])
#            
#            (data[folder])[gene[0]] = [allele_list]
#            
#    return data, numSamples, startingSamp

'''
Return data(dictionary) with all proportions rounded to 3 dp

Input: data(dictionary) loaded previously
'''
def roundProp(data):
    roundedData = dict()
    locus = list(data[data.keys()[0]].columns.values)   
    
    for sample in data:
        #sample is a key, sampleDF is a dataframe
        sampleDF = data[sample]
        roundedSampleDF = pd.DataFrame(columns=locus)   
        
        track=0     #track locus name
        for column in sampleDF.values[0]:   #column is a dictionary containing {variant : proportion}
            prop = column.values()
            keys = column.keys()
            prop = [float(p)*1000 for p in prop]    
            prop = np.array(prop)
            
            frac = prop - np.array([int(p) for p in prop])
            numFracRoundUp = int(round(sum(frac)))
            
            sortedIndices = frac.argsort()[::-1]    #descending order
            roundUpIndices = sortedIndices[0:numFracRoundUp]    #grab first numFracRoundUp indices to round up
            mask = np.zeros_like(prop, dtype=bool)  #use to identify which to round up and round down
            mask[roundUpIndices] = True
            prop[mask] = np.ceil(prop[mask])
            prop[~mask] = np.floor(prop[~mask])
            prop = prop/1000.0
            prop = ["{:.3f}".format(i) for i in prop]   #convert back to string with 3dp
            
            #reconstruct dataframe having same format as input data
            roundedDict = dict(itertools.izip(keys,prop))
            roundedSampleDF[locus[track]] = [roundedDict]
            track = track+1
            
        roundedData[sample] = roundedSampleDF
                   
    return roundedData    
            
            
'''
Check the proportions at each locus for each sample sum to 1. Raise systemExit exception and report corresponding (sample, locus) which do not sum to 1

Input: 
data: dictionary, data file loaded previously
'''            
def checkProp(data, propFormat):
    print("....Checking sum of proportions at each locus for each sample are equal to 1....")
    report = list()
    
    for sample in data:
        sampleDF = data[sample]
        
        for column in sampleDF.values[0]:
            match = re.search("(.*)_", (column.keys())[0])
            locus = match.group(1)  #get locus name
            prop = column.values()
            prop = [float(p) for p in prop]
            
            if round(sum(prop), 4) != float(propFormat):
                report.append((sample, locus))
             
    try:
        if report:
            raise SystemExit
                    
    except:
        sys.exit("The following (sample, locus) pairs have sum of proportions not equal to 1.\n {0}".format(report))
                
    print("Sum of proportions at each locus for each sample are equal to 1\n")
    
    
'''Return the unique combinations of variants at all loci across all samples
data: Data file preloaded previously

def uniqueCombinations(data):
    uniqueStrains = list()
    for sample in data:
        #key = sample name, value = dataframe
        sampleDF = data[sample]
        variantsAtAllLoci = list()
        #only one row in the dataframe
        for column in sampleDF.values[0]:
            variantsAtAllLoci.append(column.keys())
        
        combination = itertools.product(*variantsAtAllLoci)
        for strain in combination:
            uniqueStrains.append(strain)
            
    uniqueStrains = list(set(uniqueStrains))
    return uniqueStrains
'''

'''Returns a dataframe of combinations of the variants at all loci, the matrix representation of the combination
and the sample that contains it. Also, it returns the total number of combinations that each sample has

Input
data: dictionary, information of all samples preloaded previously
numLoci: number of loci
loci: list of locus
'''
def returnCombinationsAndNumComb(data, numLoci, loci):
    strains = list()
    numOfComb = dict()
    previousNum = 0
#    numAllele = dict()
    for sample in data:
        #key = sample name, value = dataframe
        sampleDF = data[sample]
        variantsAtAllLoci = list()
        
        #Use for getting matrix representation of combinations
        #only one row in the dataframe
        for locus in sampleDF.columns.tolist():
            variantsAtAllLoci.append(sampleDF[locus][0].keys())
#            numAllele[(sample, locus)] = len(sampleDF[locus][0].keys())
        
        combination = itertools.product(*variantsAtAllLoci) #produce combinations
        combinationIndices = [list(comb) for comb in itertools.product(*[range(len(var)) for var in variantsAtAllLoci])]
        for strain,strainIndex in itertools.izip(combination,combinationIndices):
            temp = list(strain)
            temp.append(sample) #add sample name
            #matrixRep = np.zeros(shape=(maxNumVar, numLoci))
            #matrixRep[strainIndex, np.arange(numLoci)] = 1
            #temp.append(matrixRep)
            strains.append(temp)
            
        #Get the total number of combinations for each sample
        numOfComb[sample] = len(strains) - previousNum
        previousNum = numOfComb[sample]
        
    strains = pd.DataFrame(strains, columns=(loci+['Sample']))
    return strains, numOfComb

'''Returns a dictionary where keys=sample name and values=numpy matrix representation of the proportions for the sample
Input
data: dictionary, information of all samples preloaded previously
numLoci: number of loci


def returnProportions(data, numLoci):
    proportions = dict()
    
    for sample in data:
        sampleDF = data[sample]
        maxNumVar = 0
        
        for column in sampleDF.values[0]:
            numVar = len(column.keys())
            if numVar > maxNumVar:
                maxNumVar = numVar
                
        proportions[sample] = np.zeros(shape=(maxNumVar, numLoci))
        
        i = 0
        for column in sampleDF.values[0]:
            numVar = len(column.keys())
            (proportions[sample])[0:numVar, i] = np.transpose(np.asarray(column.values(), dtype='float64'))
            i = i+1          
    
    return proportions            
'''

'''
Return a data frame with variant, locus, variant's proportions, the sample referring to as columns
Input:
        data: dictionary, data loaded previously
'''
def returnVarAndProportions(data):
    varAndProp = pd.DataFrame(columns=["Variant", "Locus", "Proportion", "Sample"])
    var=list()
    prop=list()
    loc = list()
    samp = list()
    
    for sample in data:
        #d is a dataframe
        d = data[sample]
        for column in d:
            #d[column][0] is a dictionary {var: prop}
            var.append((d[column])[0].keys())
            prop.append((d[column])[0].values())
            samp.append([sample]*len((d[column])[0].keys()))
            loc.append([column]*len((d[column])[0].keys()))
          
    var = [item for sublist in var for item in sublist]
    prop = [item for sublist in prop for item in sublist]
    prop = [float(i) for i in prop]
    loc = [item for sublist in loc for item in sublist]
    samp = [item for sublist in samp for item in sublist]

    varAndProp["Variant"] = var
    varAndProp["Locus"] = loc
    varAndProp["Proportion"] = prop
    varAndProp["Sample"] = samp
              
    return varAndProp

'''
Return a data frame with locus as columns, and each column contains a dictionary(wrapped by a list) which shows the indices of unique strains that each variants map to
Input:
    uniqueStr: a dataframe with unique strain types
    loci : a list containing locus name
                                                                         
def mapVarToStrain(uniqueStr, loci):
    varToStr = pd.DataFrame(columns=loci)
    for name in loci:
        uniqueVar = uniqueStr.drop_duplicates(name)
        uniqueVar = (uniqueVar[name]).reset_index(drop=True)
        var = dict()
        for row in uniqueVar:
            print row
            strList = uniqueStr[uniqueStr[name].isin([row])]["ST"]
            var[row] = strList
               
        varToStr[name] = [var]
        
    return varToStr
'''

'''
Return a data frame with locus as columns, and each column contains a dictionary(wrapped by a list) which shows the indices of unique strains(ST) that each key=(variants, sample) maps to
Input:
    strain: a dataframe with variants combination
    loci : a list containing locus name
    start: starting sample number
    numSamp: number of samples
'''                                                                           
def mapVarAndSampleToStrain(strain, loci, allSamples):
    varToStr = pd.DataFrame(columns=loci)
    for name in loci:
        varDF = strain.drop_duplicates(name)    #unique variants at a locus
        varDF = (varDF[name]).reset_index(drop=True)
        varDict = dict()
        for sample, var in list(itertools.product(allSamples, varDF.tolist())):
            #grab strain types which (sample, var) maps to
            strList = strain[(strain[name]==var) & (strain["Sample"]==sample)]["ST"]
            if len(strList) != 0:   #if not empty
                varDict[(var, sample)] = strList
               
        varToStr[name] = [varDict]
        
    return varToStr

#def writeInfoToCsv():
#    proportionWeightDecVarDF.to_csv('proportion.csv')
#    strainWeightDecVarDF.to_csv('strain.csv')
#    varAndProp.to_csv('error.csv')
#    
#    piDotCombDF = pd.DataFrame([propDecVar[0] for propDecVar in piDotComb])
#    piDotCombDF.to_csv('piDotComb.csv')
#    pd.DataFrame(propConstrRHS).to_csv('propConstrRHS.csv')


#def hammingWeightsStrain(df, loci,reference):
#    ref_copy = reference[:]
#    
#    for l in loci:
#        df.loc[:,l] = df.loc[:,l].str.split('_').str[1]
#        ref_copy.loc[:,l] = ref_copy.loc[:,l].str.split('_').str[1]
#        
#    matrix_ref = ref_copy[loci].as_matrix()
#            
#    index_hamming_dict = {ind:1.0 for ind in df.index}
#    
#    for row_weight in df[loci].itertuples():
#        ind = row_weight[0]
#        vec = list(row_weight[1:])
#        
#        for row_ref in matrix_ref:
#            hd = hamming(vec, row_ref)
#            
#            if hd < 0.5:
#                index_hamming_dict[ind] = 0.5
#                break
#            
#    return pd.DataFrame(data=index_hamming_dict.values(), index=index_hamming_dict.keys())
  

def computeDist_byGene(allele1, allele2, gene, distMat_dict):
    
    try:
        d = distMat_dict[gene][0].loc[allele1, allele2]
    except KeyError:
        print("{0} or {1} allele is not in the distance matrix dataframe".format(allele1, allele2))
        d = -1
        
    return d

'''
    Compute the minimum distance between a strain and strains in library
    Assuming strains in reference are indexed with [gene]_[index] i.e. clpA_1, clpX_5,... As we use a df.isin() method, order
    does not take into account.
'''
def computeStrDist(strain, distMat_dict, loci, reference):
    sorted_strain = sorted(strain)
    
    minDist = float('Inf')
    for row in reference[loci].itertuples(index=False):
        sorted_ref = sorted(list(row))
        assert(len(sorted_ref) == len(sorted_strain)), "Strains are of different length"
        
        temp_dist = 0
        innerLoopBreak = False
        for i in range(len(sorted_ref)):
            gene = sorted_ref[i].split("_")[0]
            d = computeDist_byGene(sorted_ref[i], sorted_strain[i], gene, distMat_dict)
            
            if d == -1:
                innerLoopBreak = True
                break
            else:
                temp_dist += d
                
        if innerLoopBreak == True:
            pass
        else:
            if temp_dist < minDist:
                minDist = temp_dist
                
    return minDist

  
'''
    Compute the weights for each strain. 0 if existing, min d(s,s') for novel strain s where s' is an existing strain
'''
def computeMinDistWeights(strainWeightDecVarDF, distMat_dict, loci, reference):
    str_df = strainWeightDecVarDF[loci]
    weight_list = list()
    
    for strain in str_df.itertuples(index=False):
        if reference[loci].isin( list(strain) ).all(axis=1).sum() == 1:
            weight_list.append(0)
        else:
            weight_list.append( computeStrDist(list(strain), distMat_dict, loci, reference) )
        
    return weight_list

'''
    Assuming each edit distance matrix is named as: editDistanceMatrix_[gene].csv
'''
def returnDistMat_dict(pathToDistMat, loci):
    distMat_dict = dict()
    
    for l in loci:
        temp_df = pd.read_csv( os.path.join(pathToDistMat, 'editDistanceMatrix_{}.csv'.format(l)), sep=",").set_index("level_0")
        distMat_dict[l] = [ temp_df ]
        

    return distMat_dict

'''
Predict strains using MILP and output in csv file
Input:
    dataPath, absolute path to directory containing samples' alleles and proportions
    pathToDistMat, path to directory containing edit distances matrix for each gene
    refStrains, path to strain_ref.txt
    outputPath, path to output csv file
    loci, list of locus
    objectiveOption, "all" means all objective components and "noPropAndErr" means omitting proportion and error terms
    globalILP_option, "all" if running on all samples
'''

def strainSolver(dataPath, refStrains, outputPath, objectiveOption, globalILP_option='all', timelimit=600, gap=8,
                 loci=["clpA","clpX","nifS","pepX","pyrG","recG","rplB","uvrA"],pathToDistMat=None):
    ''' ============================================== Data handling ====================================================== '''
    #paramaters
    propFormat = 1   #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    data, numSamples = readData(dataPath,loci, globalILP_option)
    newNameToOriName = dict()
    namingIndex=1
    for i in sorted(data.keys()):
        newNameToOriName["s{}".format(namingIndex)] = i
        data["s{}".format(namingIndex)] = data.pop(i) 
        namingIndex += 1
    reference = pd.read_csv(refStrains,sep="\t",usecols=range(1,numLoci+1))
    lociNames = list(reference.columns.values)
    numReference = reference.shape[0]
    allSamples = data.keys()
    
    #check proportions sum to 100
    checkProp(data, propFormat)
    
    #round the proportions to 3 decimal places
    data = roundProp(data)
    
    #As reference only contains numbers as entries, add gene name to the variants for better identification
    for name in lociNames:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
        
    #Get proportions of variants at different locus for each sample
    varAndProp = returnVarAndProportions(data)
    
    #Get the combinations at all loci across all samples
    strains, numOfComb = returnCombinationsAndNumComb(data, numLoci, loci)
    uniqueStrains = strains.drop_duplicates(loci)
    uniqueStrains = (uniqueStrains[loci]).reset_index(drop=True)
    uniqueStrains["ST"] = uniqueStrains.index.values + 1    #assign indices for strains or each unique combinations
    strains = strains.merge(uniqueStrains, indicator=True, how="left")    #assign the index to the combinations(as strain data frame contains duplicated rows)
    strains = strains.drop("_merge",1)
    
    #For each variants, get a mapping of which strains it maps to
    varSampToST = mapVarAndSampleToStrain(strains, loci, allSamples)
    
    #weights and decision variables for proportion of strains. weight=0 if the strain is in reference, otherwise =1. Notice there will be duplications of strain types
    #here because for proportions, we consider sample by sample rather than unique strain types
    proportionWeightDecVarDF = strains.merge(reference, indicator=True, how="left")
    proportionWeightDecVarDF["_merge"].replace(to_replace="both", value=0, inplace=True)
    proportionWeightDecVarDF["_merge"].replace(to_replace="left_only", value=1, inplace=True)
    proportionWeightDecVarDF = proportionWeightDecVarDF.rename(columns = {"_merge":"Weights"})
    
    #Add proportion decision variable names
    proportionWeightDecVarDF["Decision Variable"] = np.nan
    
    for samp in allSamples:
        thisSample = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Sample']
        propNameTemp = ["pi_%s_%d" %t for t in itertools.izip(thisSample, range(1,1+thisSample.shape[0]))]
        #shorter name as CPLEX can't hold name with >16 char. Use last 3 digits of sample name to name decision variables i.e. SRR2034333 -> use 333
        propNameTemp = [ele.replace("pi_{}".format(samp), "pi_{}".format(samp[-3:])) for ele in propNameTemp]   
        proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp, 'Decision Variable'] = propNameTemp
        
    #weights and decision variables for unique strain types, weight=0 if strain is in reference, otherwise=1. no duplications
    strainWeightDecVarDF = proportionWeightDecVarDF.drop_duplicates(loci)
    #temp_changeWeightDF = hammingWeightsStrain(strainWeightDecVarDF[strainWeightDecVarDF["Weights"] == 1], loci, reference)
    #strainWeightDecVarDF.loc[temp_changeWeightDF.index, "Weights"] = temp_changeWeightDF.values
    retainCol = loci + ['Weights', 'ST']
    strainWeightDecVarDF = strainWeightDecVarDF[retainCol].reset_index(drop=True)
    
    if pathToDistMat == None:
        pass
    else:
        distMat_dict = returnDistMat_dict(pathToDistMat, loci)
        minDist_W = computeMinDistWeights(strainWeightDecVarDF, distMat_dict, loci, reference)
        max_ed = max(minDist_W)
        minDist_W = [float(i)/max_ed for i in minDist_W]
        strainWeightDecVarDF["Weights"] = minDist_W
        
    strainWeightDecVarDF["Decision Variable"] = ["a{}".format(i) for i in range(1, strainWeightDecVarDF.shape[0] + 1)]
    
    '''==================================== Forming ILP here ================================================'''
    #Form a CPLEX model
    model = cplex.Cplex()
    
    #Some bound on cplex solver when gap finds it hard to converge
    timelim_cb = model.register_callback(TimeLimitCallback)
    timelim_cb.starttime = model.get_time()
    timelim_cb.timelimit = timelimit
    timelim_cb.acceptablegap = gap
    timelim_cb.aborted = False
    
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=strainWeightDecVarDF['Weights'].values.tolist(), names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
    #add proportions decision variables
    if objectiveOption == "noPropAndErr" or "noProp":
        model.variables.add(lb=[0]*proportionWeightDecVarDF.shape[0], ub=[propFormat]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))
    else:
        model.variables.add(obj=[i for i in proportionWeightDecVarDF['Weights'].values.tolist()],lb=[0]*proportionWeightDecVarDF.shape[0],ub=[propFormat]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))
        
    #add linear constraints such that for each sample, the sum of the proportions of its variants combination = 1
    propVarSumTo1 = list()
    
    for samp in allSamples:
        temp = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Decision Variable'].tolist()        
        propVarSumTo1.append([temp, [1]* len(temp)])
        
    model.linear_constraints.add(lin_expr=propVarSumTo1, rhs=[propFormat]*len(propVarSumTo1), senses=["E"]*len(propVarSumTo1), names=["c{0}".format(i+1) for i in range(len(propVarSumTo1))])                                                                 
    
    #add linear constraints such that each decision variable a_i must be at least pi_jk in which pi_jk is the proportion of V_jk and V_jk=a_i
    #By this, if we use any of the pi, we force a_i to be 1
    indicLargerPropDF = pd.DataFrame(columns=["ST","Indicator"])
    indicLargerPropDF["ST"] = strainWeightDecVarDF["ST"]
    indicLargerPropDF["Indicator"] = strainWeightDecVarDF["Decision Variable"]
    indicLargerPropDF = (indicLargerPropDF.merge(proportionWeightDecVarDF, indicator=True, how="left", on="ST"))[["ST","Indicator","Decision Variable"]]
    indicLargerPropDF.rename(columns={"Decision Variable": "Proportion Variable"}, inplace=True)
    indicMinusProp = list()
    for i,pi in itertools.izip(indicLargerPropDF["Indicator"].tolist(), indicLargerPropDF["Proportion Variable"].tolist()):
        indicMinusProp.append([[i, pi],[propFormat, -1]])  
    
    model.linear_constraints.add(lin_expr=indicMinusProp, rhs=[0]*len(indicMinusProp), senses=["G"]*len(indicMinusProp), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusProp))] )
    
    #Also, add linear constraints such that a_i - average of pi_jk <= 0.999. Otherwise will have case that a_i=1 and for all pi_jk, pi_jk=0
    indicMinusAvgPropLess1_DF = indicLargerPropDF.groupby("Indicator")["Proportion Variable"].apply(list).reset_index()
    indic = indicMinusAvgPropLess1_DF["Indicator"].tolist()
    pV = indicMinusAvgPropLess1_DF["Proportion Variable"].tolist()
    indicMinusAvgPropLess1_LHS = list()
    
    for i in range(len(indic)):
        a_i = indic[i]
        pi_i = pV[i]
        temp = list()
        size = len(pi_i)
        temp.append(a_i)
        coef = list()
        coef.append(propFormat)
        
        for j in range(size):
            temp.append(pi_i[j])
            coef.append(-1.0/size)
            
        indicMinusAvgPropLess1_LHS.append([temp, coef])
    
    tolerance = 0.01*propFormat*0.01     #how much tolerance we set for the upper bound    
    model.linear_constraints.add(lin_expr=indicMinusAvgPropLess1_LHS, rhs=[propFormat - tolerance]*len(indicMinusAvgPropLess1_LHS), senses=["L"]*len(indicMinusAvgPropLess1_LHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusAvgPropLess1_LHS))])
    model.linear_constraints.add(lin_expr=indicMinusAvgPropLess1_LHS, rhs=[0]*len(indicMinusAvgPropLess1_LHS), senses=["G"]*len(indicMinusAvgPropLess1_LHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusAvgPropLess1_LHS))])
    
    #add error variables and linear constraints related to error terms
    #create error variable names
    varAndProp["Decision Variable"] = ["d_{}_".format(samp[-3:]) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Decision Variable"] = varAndProp["Decision Variable"] + varAndProp["Variant"]
    
    #add error variables
    #errorThres = 0.1
    if objectiveOption == "noPropAndErr":
        model.variables.add(names=varAndProp["Decision Variable"].tolist(), lb= [0]*varAndProp.shape[0], ub= [errorThres]*varAndProp.shape[0], types=[model.variables.type.continuous]*varAndProp.shape[0])
    else:
        model.variables.add(obj=[1]*varAndProp.shape[0], names=varAndProp["Decision Variable"].tolist(), lb= [0]*varAndProp.shape[0], ub= [errorThres]*varAndProp.shape[0], types=[model.variables.type.continuous]*varAndProp.shape[0])
    
    #add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
    piDotComb = list()
    piDotComb_2 = list()
    propConstrRHS = list()
    for locusName in varSampToST:
        temp=list()
        varSampToSTDict = varSampToST[locusName][0]
        for (var, sample) in varSampToSTDict:
            strainTypes = varSampToSTDict[(var, sample)]
            propDecVar = proportionWeightDecVarDF[(proportionWeightDecVarDF["ST"].isin(strainTypes)) & (proportionWeightDecVarDF["Sample"] == "{}".format(sample))]["Decision Variable"]
            errorDecVar = varAndProp[(varAndProp["Variant"] == var) & (varAndProp["Sample"] == sample)]["Decision Variable"]
            propConstrRHS.append(  float( ( (data["{}".format(sample)])[locusName][0] )[var] )  )
            piDotComb.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [-1]])
            piDotComb_2.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [1]])
           
    model.linear_constraints.add(lin_expr=piDotComb, rhs=propConstrRHS, senses=["L"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])     
    model.linear_constraints.add(lin_expr=piDotComb_2, rhs=propConstrRHS, senses=["G"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))]) 
    
    #Export some info for MATLAB use
    #writeInfoToCsv()
    
    ''' ================================== Solve ILP ========================================== '''
    #model.write("borreliaLP.lp")
#    model.set_results_stream(None)
    #Some bound on cplex solver when gap finds it hard to converge
    timelim_cb = model.register_callback(TimeLimitCallback)
    timelim_cb.starttime = model.get_time()
    timelim_cb.timelimit = int(timelimit)
    timelim_cb.acceptablegap = float(gap)
    timelim_cb.aborted = False
    model.solve()
    
    #options for searching more optimal solutions
    #model.parameters.mip.pool.capacity.set(10)
#    model.parameters.mip.pool.intensity.set(4)
    #model.parameters.mip.limits.populate.set(2100000000)
#    model.parameters.mip.pool.absgap.set(0)
#    model.parameters.mip.pool.replace.set(1)
#    model.populate_solution_pool()
    
    objvalue = model.solution.get_objective_value()
    varNames = model.variables.get_names()
    varValues = model.solution.get_values(varNames)
    conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
    conclusion["Decision Variable"] = varNames
    conclusion["Value"] = varValues
    print conclusion
    strainInfo = conclusion.merge(strainWeightDecVarDF[strainWeightDecVarDF["Decision Variable"].isin(varNames)])
    strainInfo["New/Existing"] = ["Existing" if w==0 else "New" for w in strainInfo["Weights"].tolist()]
    strainsNeeded = (strainInfo[strainInfo["Value"] > 0.9][loci + ["ST", "New/Existing"]])
    errorVariables = conclusion[conclusion["Decision Variable"].str.startswith("d_")]
    errorVariables = errorVariables[errorVariables["Value"] != 0.0]
    strainVariables = conclusion[conclusion["Decision Variable"].str.startswith("a")]
    strainVariables = strainVariables[strainVariables["Value"] > 0.0]
    strainVariables = strainVariables.merge(strainWeightDecVarDF)
    strainVariables = strainVariables[strainVariables["Weights"] > 0]
    propVariables = conclusion[conclusion["Decision Variable"].str.startswith("pi_")]
    propVariables = propVariables[propVariables["Value"] > 0.0]
    propVariables = propVariables.merge(proportionWeightDecVarDF)
    propVariables = propVariables[propVariables["Weights"] > 0]
    print("(Strain, Proportion, Error): ({0},{1},{2})".format(strainVariables["Value"].sum(), propVariables["Value"].sum(), errorVariables["Value"].sum()))
    print("Total: {}".format(objvalue))
    
    #output indices of all strains (New/Existing)
    allStr = strainWeightDecVarDF[["ST", "Weights"] + loci]
    allStr["New/Existing"] = ["Existing" if w==0 else "New" for w in allStr["Weights"].tolist()]
    allStr.drop("Weights", 1, inplace=True)
    allStr.to_csv("{0}/indexedStrains.csv".format(outputPath))
    
#    for i in range(model.solution.pool.get_num()):
#        objvalue = model.solution.pool.get_objective_value(i)
#        varNames = model.variables.get_names()
#        varValues = model.solution.pool.get_values(i,varNames)
#        conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
#        conclusion["Decision Variable"] = varNames
#        conclusion["Value"] = varValues
#        strainInfo = conclusion.merge(strainWeightDecVarDF[strainWeightDecVarDF["Decision Variable"].isin(varNames)])
#        strainInfo["New/Existing"] = ["Existing" if w==0 else "New" for w in strainInfo["Weights"].tolist()]
#        strainsNeeded = (strainInfo[strainInfo["Value"] > 0.9][loci + ["ST", "New/Existing"]])
#        print strainsNeeded
    output_dict = dict()
    for samp in allSamples:
        output = proportionWeightDecVarDF[proportionWeightDecVarDF["Sample"] == samp].merge(strainsNeeded).drop(["Weights", "Sample"],1)
        output["Proportion"] = model.solution.get_values(output["Decision Variable"].tolist())
        output = output[output["Proportion"] > 0.0]
	output.drop("Decision Variable", axis=1, inplace=True)
        output = output[["ST", "New/Existing"]+loci+["Proportion"]]
        #print output
        temp = output
        output_dict[samp] = pd.DataFrame(temp)
        output.to_csv("{0}/{1}_strainsAndProportions.csv".format(outputPath, newNameToOriName[samp]))
   
    return output
'''
Solve the pure ILP instance of strain prediction
Input:
    dataPath, path to directory containing samples' alleles and proportions
    refStrains, path to strain_ref.txt
    loci, a list of locus
    globalILP_option, "all" if all samples
Output:
    solution_dict, dictionary where key=indices, value=dataframe related to that solution(information such as alleles at each locus)
    objective_dict, dictionary where key=indices, value=objective value of the i-th solution
    data, dictionary where key=sample name, value=dataframe which contains information about alleles and proportion at each locus
    strains, dataframe of unique strains for later use of localLP
    
'''
def minNewStrain(dataPath, refStrains, loci, globalILP_option):
    #paramaters
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    data, numSamples = readData(dataPath,loci, globalILP_option)
    newNameToOriName = dict()
    namingIndex=1
    for i in sorted(data.keys()):
        newNameToOriName["s{}".format(namingIndex)] = i
        data["s{}".format(namingIndex)] = data.pop(i) 
        namingIndex += 1
    
    allSamples = data.keys()
    reference = pd.read_csv(refStrains,sep="\t",usecols=range(1,numLoci+1))
    lociNames = list(reference.columns.values)
    
    #As reference only contains numbers as entries, add gene name to the variants for better identification
    for name in lociNames:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
    
    #Get the combinations at all loci across all samples
    strains, numOfComb = returnCombinationsAndNumComb(data, numLoci, loci)
    uniqueStrains = strains.drop_duplicates(loci)
    uniqueStrains = (uniqueStrains[loci]).reset_index(drop=True)
    uniqueStrains["ST"] = uniqueStrains.index.values + 1    #assign indices for strains or each unique combinations
    strains = strains.merge(uniqueStrains, indicator=True, how="left")    #assign the index to the combinations(as strain data frame contains duplicated rows)
    strains = strains.drop("_merge",1)
    
    #weights and decision variables for proportion of strains. weight=0 if the strain is in reference, otherwise =1. Notice there will be duplications of strain types
    #here because for proportions, we consider sample by sample rather than unique strain types
    strainWeightDecVarDF = strains.merge(reference, indicator=True, how="left")
    strainWeightDecVarDF["_merge"].replace(to_replace="both", value=0, inplace=True)
    strainWeightDecVarDF["_merge"].replace(to_replace="left_only", value=1, inplace=True)
    strainWeightDecVarDF = strainWeightDecVarDF.rename(columns = {"_merge":"Weights"})
    strainWeightDecVarDF = strainWeightDecVarDF.drop_duplicates(loci)
    retainCol = loci + ['Weights', 'ST']
    strainWeightDecVarDF = strainWeightDecVarDF[retainCol].reset_index(drop=True)
    strainWeightDecVarDF["Decision Variable"] = ["a{}".format(i) for i in range(1, strainWeightDecVarDF.shape[0] + 1)]
    
    #Relate sample and strain decision variable 
    samp_decVar_DF = strains.merge(strainWeightDecVarDF, how="left")[loci+["Sample", "Decision Variable", "ST"]]
    
    #For each allele, get a mapping of which strains it maps to
    varSampToST = mapVarAndSampleToStrain(samp_decVar_DF, loci, allSamples)
    
    '''==================================== Forming ILP here ================================================'''
    #Form a CPLEX model
    model = cplex.Cplex()
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=strainWeightDecVarDF['Weights'].values.tolist(), names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
    
    #Add linear constraints where strains chosen are able to describe all alleles seen in all samples
    descAllAlleleLHS = list()
    for locusName in varSampToST:
        varSampToSTDict = varSampToST[locusName][0]
        for (var, sample) in varSampToSTDict:
            strainTypes = varSampToSTDict[(var, sample)]
            strainDecVar = samp_decVar_DF[(samp_decVar_DF["ST"].isin(strainTypes)) & (samp_decVar_DF["Sample"] == "{}".format(sample))]["Decision Variable"].tolist()
            descAllAlleleLHS.append([strainDecVar, [1]*len(strainDecVar)])
    
    model.linear_constraints.add(lin_expr=descAllAlleleLHS, rhs=[1]*len(descAllAlleleLHS), senses=["G"]*len(descAllAlleleLHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(descAllAlleleLHS))])

#    model.solve()
#    model.write("a.lp")
#    samp_decVar_DF.to_csv("a.csv")
    
    #options for searching more optimal solutions
    #model.parameters.mip.pool.capacity.set(10)
#    model.set_results_stream(None)
    model.parameters.mip.pool.intensity.set(4)
    model.parameters.mip.limits.populate.set(50)
    model.parameters.mip.pool.absgap.set(0)
    model.parameters.mip.pool.replace.set(1)
    model.populate_solution_pool()
    
    solution_dict = dict()
    objective_dict = dict()
    for i in range(model.solution.pool.get_num()):
        objvalue = model.solution.pool.get_objective_value(i)
        objective_dict[i] = objvalue
        varNames = model.variables.get_names()
        varValues = model.solution.pool.get_values(i,varNames)
        conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
        conclusion["Decision Variable"] = varNames
        conclusion["Value"] = varValues
    
        strainInfo = conclusion.merge(strainWeightDecVarDF[strainWeightDecVarDF["Decision Variable"].isin(varNames)])
        strainInfo["New/Existing"] = ["Existing" if w==0 else "New" for w in strainInfo["Weights"].tolist()]
        strainsNeeded = (strainInfo[strainInfo["Value"] > 0.9][loci + ["ST", "Weights"]])
        strainsNeeded.reset_index(drop=True, inplace=True)
        solution_dict[i] = strainsNeeded
        
#    print("Objective value: {}".format(objective_value))
    return solution_dict, objective_dict, data, strains, newNameToOriName

'''
Solve the LP instance of strain prediction
Input:
    solution, dataframe which contains alleles at each locus for a solution
    data, see minNewStrain
    strains, see minNewStrain
    reference, dataframe of existing strains
    loci, a list of locus
Output:
    objvalue, objective value of the LP for this solution
    errObj, value of error component
    propObj, value of proportion component
    sampleAndStrainProp, a dictionary where key=sample name and value=dataframe which contains information about the strains and their proportions
    feasible, indicator whether this solution is feasible

'''
def minNewStrainProp(solution, data, strains, refStrains, loci, newNameToOriName):
    #paramaters
    propFormat = 1    #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    reference = pd.read_csv(refStrains,sep="\t",usecols=range(1,numLoci+1))
    lociNames = list(reference.columns.values)
    numReference = reference.shape[0]
    allSamples = data.keys()
    
    #check proportions sum to 100
    checkProp(data, propFormat)
    
    #round the proportions to 3 decimal places
    data = roundProp(data)
    
    #As reference only contains numbers as entries, add gene name to the variants for better identification
    for name in lociNames:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
        
    #Get proportions of variants at different locus for each sample
    varAndProp = returnVarAndProportions(data)
    
    #Add propportion variables
    proportionWeightDecVarDF = strains.merge(solution, how='left', indicator=True)
    proportionWeightDecVarDF = proportionWeightDecVarDF[proportionWeightDecVarDF["_merge"] == "both"]
    proportionWeightDecVarDF.drop(["_merge"], axis=1, inplace=True)
    proportionWeightDecVarDF.reset_index(drop=True, inplace=True)
    
    #For each variants, get a mapping of which strains it maps to. Only consider those strains in given solution
    varSampToST = mapVarAndSampleToStrain(proportionWeightDecVarDF[loci+["Sample", "ST"]], loci, allSamples)
    
    #Add proportion variables names
    for samp in allSamples:
        thisSample = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Sample']
        propNameTemp = ["pi_%s_%d" %t for t in itertools.izip(thisSample, range(1,1+thisSample.shape[0]))]
        #shorter name as CPLEX can't hold name with >16 char. Use last 3 digits of sample name to name decision variables i.e. SRR2034333 -> use 333
        propNameTemp = [ele.replace("pi_{}".format(samp), "pi_{}".format(samp)) for ele in propNameTemp]   
        proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp, 'Decision Variable'] = propNameTemp
        
    ''' ===================================== Forming LP here =================================================== '''
    #Form a CPLEX model
    model = cplex.Cplex()
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=proportionWeightDecVarDF['Weights'].values.tolist(), lb=[0]*proportionWeightDecVarDF.shape[0], ub=[propFormat]*proportionWeightDecVarDF.shape[0], names=proportionWeightDecVarDF['Decision Variable'], types = [model.variables.type.continuous]* len(proportionWeightDecVarDF['Weights'].values.tolist()))
        
    #add linear constraints such that for each sample, the sum of the proportions of its variants combination = 1
    propVarSumTo1 = list()
    
    for samp in allSamples:
        temp = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Decision Variable'].tolist()        
        propVarSumTo1.append([temp, [1]* len(temp)])
        
    model.linear_constraints.add(lin_expr=propVarSumTo1, rhs=[propFormat]*len(propVarSumTo1), senses=["E"]*len(propVarSumTo1), names=["c{0}".format(i+1) for i in range(len(propVarSumTo1))])
    
    #add error variables and linear constraints related to error terms
    #create error variable names
    varAndProp["Decision Variable"] = ["d_{}_".format(samp) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Decision Variable"] = varAndProp["Decision Variable"] + varAndProp["Variant"]
    
     #create artificial variable to minimize absolute value of error
#    varAndProp["Artificial"] = ["f_s{}_".format(samp[-3:]) for samp in varAndProp["Sample"].tolist() ]
#    varAndProp["Artificial"] = varAndProp["Artificial"] + varAndProp["Variant"]
    
    #add error variables
    #errorThres = 0.1
    model.variables.add(obj=[1]*varAndProp.shape[0], names=varAndProp["Decision Variable"].tolist(), lb=[0]*varAndProp.shape[0], ub= [errorThres]*varAndProp.shape[0], types=[model.variables.type.continuous]*varAndProp.shape[0])
#    model.variables.add(obj=[1]*varAndProp.shape[0], names=varAndProp["Artificial"].tolist(), lb=[0]*varAndProp.shape[0], ub= [0.2]*varAndProp.shape[0], types=[model.variables.type.continuous]*varAndProp.shape[0])
#    artificial_constr1 = [[[artif, err],[1,1]] for artif, err in itertools.izip(varAndProp["Artificial"].tolist(), varAndProp["Decision Variable"].tolist())]
#    artificial_constr2 = [[[artif, err],[1,-1]] for artif, err in itertools.izip(varAndProp["Artificial"].tolist(), varAndProp["Decision Variable"].tolist())]
#    model.linear_constraints.add(lin_expr=artificial_constr1, rhs=[0]*len(artificial_constr1), senses=["G"]*len(artificial_constr1), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(artificial_constr1))])
#    model.linear_constraints.add(lin_expr=artificial_constr2, rhs=[0]*len(artificial_constr2), senses=["G"]*len(artificial_constr2), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(artificial_constr2))])
    
   #add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
    piDotComb = list()
    piDotComb_2 = list()
    propConstrRHS = list()
    for locusName in varSampToST:
        temp=list()
        varSampToSTDict = varSampToST[locusName][0]
        for (var, sample) in varSampToSTDict:
            strainTypes = varSampToSTDict[(var, sample)]
            propDecVar = proportionWeightDecVarDF[(proportionWeightDecVarDF["ST"].isin(strainTypes)) & (proportionWeightDecVarDF["Sample"] == "{}".format(sample))]["Decision Variable"]
            errorDecVar = varAndProp[(varAndProp["Variant"] == var) & (varAndProp["Sample"] == sample)]["Decision Variable"]
            propConstrRHS.append(  float( ( (data["{}".format(sample)])[locusName][0] )[var] )  )
            piDotComb.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [-1]])
            piDotComb_2.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [1]])
           
    model.linear_constraints.add(lin_expr=piDotComb, rhs=propConstrRHS, senses=["L"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])           
    model.linear_constraints.add(lin_expr=piDotComb_2, rhs=propConstrRHS, senses=["G"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])                                                                   
    
    #error must sum to 0
#    errorSumTo0 = list()
#    for samp, loc in list(set(itertools.izip(varAndProp["Sample"].tolist(), varAndProp["Locus"].tolist()))):
#        temp = (varAndProp[(varAndProp["Sample"] == samp) & (varAndProp["Locus"] == loc)])["Decision Variable"].tolist()
#        errorSumTo0.append([temp, [1]*len(temp)])
#        
#    model.linear_constraints.add(lin_expr=errorSumTo0, rhs=[0]*len(errorSumTo0), senses=["E"]*len(errorSumTo0), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errorSumTo0))])
#    model.variables.set_upper_bounds([(i, 0.2) for i in varAndProp["Artificial"].tolist()])
    
    ''' ==== Solve ==== '''
    model.set_problem_type(0)   #set to LP problem
#    model.set_results_stream(None)
#    model.set_error_stream(None)
#    model.write("a.lp")
    model.solve()
#    print model.solution.get_status_string()
    
    feasible = False
    if model.solution.get_status() == 1:
        objvalue = model.solution.get_objective_value()
        varNames = model.variables.get_names()
        varValues = model.solution.get_values(varNames)
        conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
        conclusion["Decision Variable"] = varNames
        conclusion["Value"] = varValues
        conclusion = conclusion[conclusion["Value"] != 0.0]
        print conclusion
        print conclusion[conclusion["Decision Variable"].str.startswith("d_")]["Value"].sum()
        errObj = conclusion[conclusion["Decision Variable"].str.startswith("d_")]["Value"].sum()
        propVariables = conclusion[conclusion["Decision Variable"].str.startswith("pi_")]
        propVariables = propVariables.merge(proportionWeightDecVarDF)
        propVariables = propVariables[propVariables["Weights"] != 0]
        propObj = propVariables["Value"].sum()
        
        sampleAndStrainProp = dict()
        for samp in allSamples:
            output = proportionWeightDecVarDF[proportionWeightDecVarDF["Sample"] == samp].merge(conclusion)
            output["New/Existing"] = ["Existing" if w==0 else "New" for w in output["Weights"].tolist()]
            output.rename(columns={"Value":"Proportion"}, inplace=True)
            output.drop("Decision Variable", axis=1, inplace=True)
            output = output[["ST", "New/Existing"]+loci+["Proportion"]]
            sampleAndStrainProp[newNameToOriName[samp]] = output
    #        print(output)
    #        output.to_csv("{0}/{1}_strainsAndProportions.csv".format(outputPath, samp))
        feasible = True
    else:
        objvalue= -1
        errObj = -1
        propObj = -1
        sampleAndStrainProp = list()

    return objvalue, errObj, propObj, sampleAndStrainProp, feasible

