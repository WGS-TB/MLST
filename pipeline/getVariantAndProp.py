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
    minVar_Qscores = list()
    for i in range(len(all_solutions)):
        Qscores.append(compute_QSum(Qmatrix.loc[reads_cov,all_solutions[i]]))
        
    print("Quality scores for all solutions: {}".format(Qscores))
    print("Solutions: {}".format(all_solutions))
    
    #The required alleles
    dataMatrix_pred = dataMatrixDF.loc[reads_cov,var_predicted]
#    minVar_solutions = [sol for sol in all_solutions if len(sol) == min(map(len,all_solutions))]
    minVar_solutions = all_solutions
    
    #Calculate quality scores for minimum alleles
    for i in range(len(minVar_solutions)):
        minVar_Qscores.append(compute_QSum(Qmatrix.loc[reads_cov, minVar_solutions[i]]))
        
    print("Quality scores for solutions with minimum alleles: {}".format(minVar_Qscores))
    print("Minimum alleles solutions: {}".format(minVar_solutions))
    
#    Return solutions with minimum quality scores
    min_qscore = min(Qscores)
    minQscore_sol = [minVar_solutions[i] for i in range(len(minVar_solutions)) if Qscores[i] == min_qscore]
    print("Min quality score: {}".format(min_qscore))
    print("Solutions with minimum quality score: {}".format(minQscore_sol))
    if len(minQscore_sol) > 1:
        print("@@@@@@@@@@@@ More than 1 solution having minimum quality score @@@@@@@@@@@@")
    
    #Likelihood approach 
     
    #compute proportions
#    prop = compute_proportions(dataMatrix_pred)
#    pred_prop = create_dictionary(var_predicted, prop)
    
    #score list and proportions
#    score_list = list()
#    min_score = sys.maxint
#    
#    for i in range(len(minVar_solutions)):
#        score = compute_likelihood(dataMatrixDF.loc[reads_cov, minVar_solutions[i]],6)
#        score_list.append(score)
#        
#        if score <= min_score:
#            min_score = score
#         
#
#    minLikeli = np.argmin(score_list)
#    var_predicted = all_solutions[minLikeli]
#    minLike_sol = [minVar_solutions[i] for i in range(len(minVar_solutions)) if score_list[i] == min_score]
#    minIndices = np.where(np.array(score_list) == np.array(score_list).min())
#    print("Likelihood score for all solutions: {}".format(score_list))
#    print("Solutions: {}".format(minVar_solutions))
#    print("Minimum neg log likelihood: {}".format(min_score))
#    print("Solutions having minimum neg log likelihood: {}".format(minLike_sol))
#    
#    if len(minIndices) > 1:
#        print("More than 1 solution having min likelihood")
    #Give some names to the solutions for further identifications and get the indices for sorted likelihood list
    #Sort quality score and produce some graphs
#    sortedIndex_score_list = np.argsort(score_list)
#    sortedIndex_qscore_minVarSol = np.argsort(Qscores)
#    likelihood_score_dict = dict()
#    qscores_minVarSol_dict = dict()
#    sol_name_dict = dict()
    
#    for i in range(len(minVar_solutions)):
#        sol_name_dict["sol_{}".format(i)] = minVar_solutions[sortedIndex_score_list[i]]
#        likelihood_score_dict["sol_{}".format(i)] = score_list[sortedIndex_score_list[i]]
#        qscores_minVarSol_dict["sol_{}".format(i)] = Qscores[sortedIndex_qscore_minVarSol[i]]
            
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

def maxExistingStr(sample, loci, gene_solProp_dict, reference, objectiveOption):
    genesDF = pd.DataFrame(columns=loci)
    
    for gene in loci:
        genesDF[gene] = [gene_solProp_dict[gene]]
    
    data = dict()
    data[sample] = genesDF
    data = gsp.roundProp(data)
    
    ''' ============================================== Data handling ====================================================== '''
    #paramaters
    propFormat = 1    #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    allSamples = data.keys()
        
    #Get proportions of variants at different locus for each sample
    varAndProp = gsp.returnVarAndProportions(data)
    
    #Get the combinations at all loci across all samples
    strains, numOfComb, maxAllele_dict = gsp.returnCombinationsAndNumComb(data, numLoci, loci)
#    numOfComb = strainAndNumComb[1]
    uniqueStrains = strains.drop_duplicates(loci)
    uniqueStrains = (uniqueStrains[loci]).reset_index(drop=True)
    uniqueStrains["ST"] = uniqueStrains.index.values + 1    #assign indices for strains or each unique combinations
    strains = strains.merge(uniqueStrains, indicator=True, how="left")    #assign the index to the combinations(as strain data frame contains duplicated rows)
    strains = strains.drop("_merge",1)
    
    #For each variants, get a mapping of which strains it maps to
    varSampToST = gsp.mapVarAndSampleToStrain(strains, loci, allSamples)
    
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
        propNameTemp = [ele.replace("pi_{}".format(samp), "pi_s{}".format(samp[-3:])) for ele in propNameTemp]   
        proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp, 'Decision Variable'] = propNameTemp
        
    #weights and decision variables for unique strain types, weight=0 if strain is in reference, otherwise=1. no duplications
    strainWeightDecVarDF = proportionWeightDecVarDF.drop_duplicates(loci)
    retainCol = loci + ['Weights', 'ST']
    strainWeightDecVarDF = strainWeightDecVarDF[retainCol].reset_index(drop=True)
    strainWeightDecVarDF["Decision Variable"] = ["a{}".format(i) for i in range(1, strainWeightDecVarDF.shape[0] + 1)]
    
    '''==================================== Forming ILP here ================================================'''
    #Form a CPLEX model
    model = cplex.Cplex()
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=[i for i in strainWeightDecVarDF['Weights'].values.tolist()], names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
    #add proportions decision variables
    if objectiveOption == "noPropAndErr":
        model.variables.add(ub=[propFormat]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))
    else:
        model.variables.add(obj=proportionWeightDecVarDF['Weights'].values.tolist(),ub=[propFormat]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))
    
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
    varAndProp["Decision Variable"] = ["d_s{}_".format(samp[-3:]) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Decision Variable"] = varAndProp["Decision Variable"] + varAndProp["Variant"]
    
    #create artificial variable to minimize absolute value of error
    varAndProp["Artificial"] = ["f_s{}_".format(samp[-3:]) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Artificial"] = varAndProp["Artificial"] + varAndProp["Variant"]
    
    #add error variables
    errorThres = 1
    model.variables.add(lb=[-1*errorThres*i for i in varAndProp["Proportion"].tolist()], ub=[(1-i)*errorThres for i in varAndProp["Proportion"].tolist()], names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
    model.variables.add(obj=[1]*varAndProp.shape[0], names=varAndProp["Artificial"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
    artificial_constr1 = [[[artif, err],[1,1]] for artif, err in itertools.izip(varAndProp["Artificial"].tolist(), varAndProp["Decision Variable"].tolist())]
    artificial_constr2 = [[[artif, err],[1,-1]] for artif, err in itertools.izip(varAndProp["Artificial"].tolist(), varAndProp["Decision Variable"].tolist())]
    model.linear_constraints.add(lin_expr=artificial_constr1, rhs=[0]*len(artificial_constr1), senses=["G"]*len(artificial_constr1), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(artificial_constr1))])
    model.linear_constraints.add(lin_expr=artificial_constr2, rhs=[0]*len(artificial_constr2), senses=["G"]*len(artificial_constr2), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(artificial_constr2))])
    
    #add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
    piDotComb = list()
#    piDotComb_2 = list()
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
#            piDotComb_2.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [1]])
           
    model.linear_constraints.add(lin_expr=piDotComb, rhs=propConstrRHS, senses=["E"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])     
#    model.linear_constraints.add(lin_expr=piDotComb_2, rhs=propConstrRHS, senses=["G"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))]) 
    
    #add error variable
    #model.variables.add(lb=(-1*varAndProp["Proportion"]).tolist(), ub=(1-varAndProp["Proportion"]).tolist(), names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
#    if objectiveOption == "noPropAndErr" or objectiveOption == "noErr":
#        model.variables.add(lb=(-1*varAndProp["Proportion"]).tolist(), ub=(1-varAndProp["Proportion"]).tolist(), names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
#    else:
#        model.variables.add(obj=[1]*varAndProp.shape[0], lb=(-1*varAndProp["Proportion"]).tolist(), ub=(propFormat-varAndProp["Proportion"]).tolist(), names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
        

    #add the constraints whereby for each sample, at each locus, the sum of the error of all variants=0
    errorSumTo0 = list()
    for samp, loc in list(set(itertools.izip(varAndProp["Sample"].tolist(), varAndProp["Locus"].tolist()))):
        temp = (varAndProp[(varAndProp["Sample"] == samp) & (varAndProp["Locus"] == loc)])["Decision Variable"].tolist()
        errorSumTo0.append([temp, [1]*len(temp)])
        
    model.linear_constraints.add(lin_expr=errorSumTo0, rhs=[0]*len(errorSumTo0), senses=["E"]*len(errorSumTo0), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errorSumTo0))])
    '''
    #add the constraints which bound the error terms
    errLessSumMinProp = list()
    errLessSumMinPropRHS = list()
    errLessPropMinSum = list()
    errLessPropMinSumRHS = list()
    
    for index, row in varAndProp.iterrows():
        samp = row["Sample"]
        var = row["Variant"]
        loc = row["Locus"]
        
        err = row["Decision Variable"]
        pi = proportionWeightDecVarDF[( proportionWeightDecVarDF["Sample"]  == samp ) & (proportionWeightDecVarDF[loc] == var )]["Decision Variable"].tolist()
        prop = row["Proportion"]
        
        errLessSumMinProp.append( [[err] + pi, [i for i in itertools.chain([1],[-1]*len(pi))]] )
        errLessSumMinPropRHS.append(-1*prop)
        
        errLessPropMinSum.append( [[err] + pi, [i for i in itertools.chain([1],[1]*len(pi))]] )
        errLessPropMinSumRHS.append(prop)
    
    model.linear_constraints.add(lin_expr=errLessSumMinProp, rhs=errLessSumMinPropRHS, senses=["G"]*len(errLessSumMinProp), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errLessSumMinProp))])  
    model.linear_constraints.add(lin_expr=errLessPropMinSum, rhs=errLessPropMinSumRHS, senses=["G"]*len(errLessPropMinSum), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errLessPropMinSum))])
    '''
    #Add a known optimal objective value as constraint
    #model.linear_constraints.add(lin_expr=[ [model.variables.get_names(), [1]*len(model.variables.get_names())] ], rhs=[10], senses=["L"])
    
    #Export some info for MATLAB use
    #writeInfoToCsv()
    
    ''' ================================== Solve ILP ========================================== '''
    #model.write("borreliaLP.lp")
#    model.set_results_stream(None)
    model.solve()
#    model.write("{}.lp".format(sample))
    
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

def localILP(sample, loci, gene_solProp_dict, reference):
    genesDF = pd.DataFrame(columns=loci)
    
    for gene in loci:
        genesDF[gene] = [gene_solProp_dict[gene]]
    
    data = dict()
    data[sample] = genesDF
    data = gsp.roundProp(data)
    
     #paramaters
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #Get the combinations at all loci across all samples
    strains, numOfComb, maxAllele_dict = gsp.returnCombinationsAndNumComb(data, numLoci, loci)
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
    samp_decVar_DF = strains.merge(strainWeightDecVarDF, how="left")[["Sample", "Decision Variable"]]
    samp_decVar_dict = {samp: group["Decision Variable"].tolist() for samp, group in samp_decVar_DF.groupby("Sample")}
    
    '''==================================== Forming ILP here ================================================'''
    #Form a CPLEX model
    model = cplex.Cplex()
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=strainWeightDecVarDF['Weights'].values.tolist(), names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
    
    #Add linear constraints where strains chosen are able to describe all alleles seen in all samples
    descAllAlleleLHS = list()
    descAllAlleleRHS = list()
    for samp in maxAllele_dict.keys():
        descAllAlleleLHS.append([samp_decVar_dict[samp], [1]*len(samp_decVar_dict[samp])])
        descAllAlleleRHS.append(maxAllele_dict[samp])
        
    model.linear_constraints.add(lin_expr=descAllAlleleLHS, rhs=descAllAlleleRHS, senses=["G"]*len(descAllAlleleRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(descAllAlleleRHS))])
    
#    model.solve()
    
    #options for searching more optimal solutions
    #model.parameters.mip.pool.capacity.set(10)
    model.set_results_stream(None)
    model.parameters.mip.pool.intensity.set(4)
    model.parameters.mip.limits.populate.set(100)
    model.parameters.mip.pool.absgap.set(5)
    model.parameters.mip.pool.replace.set(2)
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
        
    return solution_dict, objective_dict, data, strains

def localLP(solution, data, strains, reference, loci):
    #paramaters
    propFormat = 1    #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    lociNames = list(reference.columns.values)
    numReference = reference.shape[0]
    allSamples = data.keys()
        
    #Get proportions of variants at different locus for each sample
    varAndProp = gsp.returnVarAndProportions(data)
    
    #Add propportion variables
    proportionWeightDecVarDF = strains.merge(solution, how='left', indicator=True)
    proportionWeightDecVarDF = proportionWeightDecVarDF[proportionWeightDecVarDF["_merge"] == "both"]
    proportionWeightDecVarDF.drop(["_merge"], axis=1, inplace=True)
    proportionWeightDecVarDF.reset_index(drop=True, inplace=True)
    
    #For each variants, get a mapping of which strains it maps to. Only consider those strains in given solution
    varSampToST = gsp.mapVarAndSampleToStrain(proportionWeightDecVarDF[loci+["Sample", "ST"]], loci, allSamples)
    
    #Add proportion variables names
    for samp in allSamples:
        thisSample = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Sample']
        propNameTemp = ["pi_%s_%d" %t for t in itertools.izip(thisSample, range(1,1+thisSample.shape[0]))]
        #shorter name as CPLEX can't hold name with >16 char. Use last 3 digits of sample name to name decision variables i.e. SRR2034333 -> use 333
        propNameTemp = [ele.replace("pi_{}".format(samp), "pi_s{}".format(samp[-3:])) for ele in propNameTemp]   
        proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp, 'Decision Variable'] = propNameTemp
        
    ''' ===================================== Forming LP here =================================================== '''
    #Form a CPLEX model
    model = cplex.Cplex()
    #minimize problem
    model.objective.set_sense(model.objective.sense.minimize)
    #add the decision variables for unqiue strain types
    model.variables.add(obj=proportionWeightDecVarDF['Weights'].values.tolist(), names=proportionWeightDecVarDF['Decision Variable'], types = [model.variables.type.continuous]* len(proportionWeightDecVarDF['Weights'].values.tolist()))
        
    #add linear constraints such that for each sample, the sum of the proportions of its variants combination = 1
    propVarSumTo1 = list()
    
    for samp in allSamples:
        temp = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Decision Variable'].tolist()        
        propVarSumTo1.append([temp, [1]* len(temp)])
        
    model.linear_constraints.add(lin_expr=propVarSumTo1, rhs=[propFormat]*len(propVarSumTo1), senses=["E"]*len(propVarSumTo1), names=["c{0}".format(i+1) for i in range(len(propVarSumTo1))])
    
    #add error variables and linear constraints related to error terms
    #create error variable names
    varAndProp["Decision Variable"] = ["d_s{}_".format(samp[-3:]) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Decision Variable"] = varAndProp["Decision Variable"] + varAndProp["Variant"]
    
     #create artificial variable to minimize absolute value of error
    varAndProp["Artificial"] = ["f_s{}_".format(samp[-3:]) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Artificial"] = varAndProp["Artificial"] + varAndProp["Variant"]
    
    #add error variables
    errorThres = 1
    model.variables.add(lb=[-1*errorThres*i for i in varAndProp["Proportion"].tolist()], ub=[(1-i)*errorThres for i in varAndProp["Proportion"].tolist()], names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
#    model.variables.add(names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
    model.variables.add(obj=[1]*varAndProp.shape[0], names=varAndProp["Artificial"].tolist(), ub= [0.15]*varAndProp.shape[0], types=[model.variables.type.continuous]*varAndProp.shape[0])
    artificial_constr1 = [[[artif, err],[1,1]] for artif, err in itertools.izip(varAndProp["Artificial"].tolist(), varAndProp["Decision Variable"].tolist())]
    artificial_constr2 = [[[artif, err],[1,-1]] for artif, err in itertools.izip(varAndProp["Artificial"].tolist(), varAndProp["Decision Variable"].tolist())]
    model.linear_constraints.add(lin_expr=artificial_constr1, rhs=[0]*len(artificial_constr1), senses=["G"]*len(artificial_constr1), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(artificial_constr1))])
    model.linear_constraints.add(lin_expr=artificial_constr2, rhs=[0]*len(artificial_constr2), senses=["G"]*len(artificial_constr2), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(artificial_constr2))])
    
   #add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
    piDotComb = list()
#    piDotComb_2 = list()
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
#            piDotComb_2.append([propDecVar.tolist() + errorDecVar.tolist(), [1]*len(propDecVar) + [1]])
           
    model.linear_constraints.add(lin_expr=piDotComb, rhs=propConstrRHS, senses=["E"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])                                                                              
    
    ''' ==== Solve ==== '''
    model.set_problem_type(0)   #set to LP problem
    model.set_results_stream(None)
    model.set_error_stream(None)
    model.solve()
    
    objvalue = model.solution.get_objective_value()

    return objvalue
    
