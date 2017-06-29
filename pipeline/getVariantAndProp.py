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
    prob_list = [0.0]*dataframe.shape[1]
    
    for row in dataframe.itertuples(index=False):
        mmInfo = [i for i in list(row) if i>=0]
        min_mm = min(mmInfo)
        numOfVar_minMm = len([i for i in list(row) if i== min_mm])
        
        for i in range(len(list(row))):
            if list(row)[i] == min_mm:
                prob_list[i] += 1/numOfVar_minMm
                
    normalize_term = 1.0/(sum(prob_list))
    prob_list = [100.0*normalize_term * i for i in prob_list]
    return prob_list 

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

#Only return solutions with minimum number of variants
def getVarAndProp(gene, tablePath, samp):
    #generate matrix
    dataMatrixDF = generate_matrix(tablePath)
    dataMatrixDF.rename(columns={'Unnamed: 0': 'Read'}, inplace=True)
    #predict variants
    pred_object_val,var_predicted,reads_cov, all_solutions, all_objective = varSolver.solver(dataMatrixDF)
    #print var_predicted
    #print all_solutions
    dataMatrix_pred = dataMatrixDF.loc[reads_cov,var_predicted]
    minVar_solutions = [sol for sol in all_solutions if len(sol) == min(map(len,all_solutions))]
    
    ''' Likelihood '''
    #compute proportions
    prop = compute_proportions(dataMatrix_pred)
    pred_prop = create_dictionary(var_predicted, prop)
    
    #score list and proportions
    score_list = list()
    min_score = sys.maxint
    min_sol_list = list()
    
    for i in range(len(minVar_solutions)):
        score = compute_likelihood(dataMatrixDF.loc[reads_cov, minVar_solutions[i]])
        score_list.append(score)
        
        if score <= min_score:
            min_score = score
            
    #Give some names to the solutions for further identifications and get the indices for sorted likelihood list
    sortedIndex_score_list = np.argsort(score_list)
    likelihood_score_dict = dict()
    sol_name_dict = dict()
    
    for i in range(len(minVar_solutions)):
        sol_name_dict["sol_{}".format(i)] = minVar_solutions[sortedIndex_score_list[i]]
        likelihood_score_dict["sol_{}".format(i)] = score_list[sortedIndex_score_list[i]]
            
    plt.figure()
    sorted_solution_namelist = ["sol_{}".format(i) for i in range(len(minVar_solutions))]
    plt.xticks(range(len(minVar_solutions)), sorted_solution_namelist, rotation=20)
    plt.scatter(range(len(minVar_solutions)), [likelihood_score_dict[name] for name in sorted_solution_namelist], s=50)
    plt.xlabel('Solution i ')
    plt.ylabel('Negative log likelihood')
    plt.savefig("{0}_{1}_sol_likelihood".format(samp, gene))
    
    minVar_minNegLog_solutions = [minVar_solutions[i] for i in range(len(minVar_solutions)) if min_score <= score_list[i] <= 1.01*min_score]
    
    ''' ====== '''
    
    #compute proportions
    #solutionsAndProp_dict is a dictionary in which the keys are just indices and values are dictionaries, with variant as key and proportion as value
    solutionsAndProp_dict = dict()
    track=0
    for sol in minVar_minNegLog_solutions:
        dataMatrix_pred = dataMatrixDF.loc[reads_cov, sol]
        prop = compute_proportions(dataMatrix_pred)
        pred_prop = create_dictionary(sol, prop)
        solutionsAndProp_dict[track] = pred_prop
        track += 1
    
    #write proportions to file
#    w = csv.writer(open(gene+'_proportions.csv', "w"))
#    for key, val in pred_prop.items():
#        w.writerow([key, val])
        
    plt.close('all')        
    return solutionsAndProp_dict

def maxExistingStr(sample, loci, gene_solProp_dict, reference):
    genesDF = pd.DataFrame(columns=loci)
    
    for gene in loci:
        genesDF[gene] = [gene_solProp_dict[gene]]
    
    data = dict()
    data[sample] = genesDF
#    data = gsp.roundProp(data)
    
    ''' ============================================== Data handling ====================================================== '''
    #paramaters
    propFormat = 100    #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    allSamples = data.keys()
        
    #Get proportions of variants at different locus for each sample
    varAndProp = gsp.returnVarAndProportions(data)
    
    #Get the combinations at all loci across all samples
    strainAndNumComb = gsp.returnCombinationsAndNumComb(data, numLoci, loci)
    strains = strainAndNumComb[0]
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
    model.variables.add(obj=strainWeightDecVarDF['Weights'].values.tolist(), names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
    #add proportions decision variables
    model.variables.add(obj=proportionWeightDecVarDF['Weights'].values.tolist(),ub=[propFormat]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))
    
    #add linear constraints such that for each sample, the sum of the proportions of its variants combination = 1
    propVarSumTo1 = list()
    
    for samp in allSamples:
        temp = (proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == samp])['Decision Variable'].tolist()        
        propVarSumTo1.append([temp, [1]* len(temp)])
        
    model.linear_constraints.add(lin_expr=propVarSumTo1, rhs=[propFormat]*len(propVarSumTo1), senses=["E"]*len(propVarSumTo1), names=["c{0}".format(i+1) for i in range(len(propVarSumTo1))])
    
    #add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
    piDotComb = list()
    propConstrRHS = list()
    for locusName in varSampToST:
        temp=list()
        varSampToSTDict = varSampToST[locusName][0]
        for (var, sample) in varSampToSTDict:
            strainTypes = varSampToSTDict[(var, sample)]
            propDecVar = proportionWeightDecVarDF[(proportionWeightDecVarDF["ST"].isin(strainTypes)) & (proportionWeightDecVarDF["Sample"] == "{}".format(sample))]["Decision Variable"]
            propConstrRHS.append(  float( ( (data["{}".format(sample)])[locusName][0] )[var] )  )
            piDotComb.append([propDecVar.tolist(), [1]*len(propDecVar)])
           
    model.linear_constraints.add(lin_expr=piDotComb, rhs=propConstrRHS, senses=["E"]*len(propConstrRHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(propConstrRHS))])                                                                         
    
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
    
    tolerance = 0.01     #how much tolerance we set for the upper bound    
    model.linear_constraints.add(lin_expr=indicMinusAvgPropLess1_LHS, rhs=[propFormat - tolerance]*len(indicMinusAvgPropLess1_LHS), senses=["L"]*len(indicMinusAvgPropLess1_LHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusAvgPropLess1_LHS))])
    model.linear_constraints.add(lin_expr=indicMinusAvgPropLess1_LHS, rhs=[0]*len(indicMinusAvgPropLess1_LHS), senses=["G"]*len(indicMinusAvgPropLess1_LHS), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusAvgPropLess1_LHS))])
    
    #add error variables and linear constraints related to error terms
    #create error variable names
    varAndProp["Decision Variable"] = ["d_s{}_".format(samp[-3:]) for samp in varAndProp["Sample"].tolist() ]
    varAndProp["Decision Variable"] = varAndProp["Decision Variable"] + varAndProp["Variant"]
              
    #add error variable
    #model.variables.add(lb=(-1*varAndProp["Proportion"]).tolist(), ub=(1-varAndProp["Proportion"]).tolist(), names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
    model.variables.add(obj=[1]*varAndProp.shape[0], lb=(-1*varAndProp["Proportion"]).tolist(), ub=(propFormat-varAndProp["Proportion"]).tolist(), names=varAndProp["Decision Variable"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])
    
    #add the constraints whereby for each sample, at each locus, the sum of the error of all variants=0
    errorSumTo0 = list()
    for samp, loc in list(set(itertools.izip(varAndProp["Sample"].tolist(), varAndProp["Locus"].tolist()))):
        temp = (varAndProp[(varAndProp["Sample"] == samp) & (varAndProp["Locus"] == loc)])["Decision Variable"].tolist()
        errorSumTo0.append([temp, [1]*len(temp)])
        
    model.linear_constraints.add(lin_expr=errorSumTo0, rhs=[0]*len(errorSumTo0), senses=["E"]*len(errorSumTo0), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errorSumTo0))])
    
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
    
    model.linear_constraints.add(lin_expr=errLessSumMinProp, rhs=errLessSumMinPropRHS, senses=["L"]*len(errLessSumMinProp), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errLessSumMinProp))])  
    model.linear_constraints.add(lin_expr=errLessPropMinSum, rhs=errLessPropMinSumRHS, senses=["L"]*len(errLessPropMinSum), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errLessPropMinSum))])
    
    #Add a known optimal objective value as constraint
    #model.linear_constraints.add(lin_expr=[ [model.variables.get_names(), [1]*len(model.variables.get_names())] ], rhs=[10], senses=["L"])
    
    #Export some info for MATLAB use
    #writeInfoToCsv()
    
    ''' ================================== Solve ILP ========================================== '''
    #model.write("borreliaLP.lp")
    model.set_results_stream(None)
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
    error = conclusion[conclusion["Decision Variable"].str.contains("d_")]
    nonZero_error = error[(error["Value"] <= -0.001) | (error["Value"] >= 0.001)]
    if nonZero_error.shape[0] != 0:
        print("^^^^^^^^^^^^^^^^^^^^^^^^^^^ error ^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        print nonZero_error
        
    objStr = conclusion[conclusion["Decision Variable"].str.contains("^a")]["Decision Variable"].tolist()
    objProp = conclusion[conclusion["Decision Variable"].str.contains("^pi")]["Decision Variable"].tolist()
    objErr = conclusion[conclusion["Decision Variable"].str.contains("^d")]["Decision Variable"].tolist()
    
    objStr_coeff = model.objective.get_linear(objStr)
    objProp_coeff = model.objective.get_linear(objProp)
    objErr_coeff = model.objective.get_linear(objErr)
    
    sum_str = sum([val*coeff for val, coeff in itertools.izip(model.solution.get_values(objStr), objStr_coeff)])
    sum_prop = sum([val*coeff for val, coeff in itertools.izip(model.solution.get_values(objProp), objProp_coeff)] )
    sum_err = sum([val*coeff for val, coeff in itertools.izip(model.solution.get_values(objErr), objErr_coeff)] )
    print("Objective value: {}".format(objvalue))
    print("Strain componenet: {}".format(sum_str))
    print("Prop componenet: {}".format(sum_prop))
    print("Error componenet: {}".format(sum_err))
    print("Sum :{}".format(sum_str+sum_prop+sum_err))
    
    return sum_str, sum_prop, sum_err
    
    
