#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:33:50 2017

@author: stanleygan

Project: Illuminating the diversity of pathogenic bacteria Borrelia Burgdorferi in tick samples

This program reads the data, construct a mixed integer linear program and solves it. The data consists
of different samples, each has 8 locus, each locus has different variants and their proportions. The 
Mixed ILP is implemented using CPLEX Python API. 

In the end, this program will output the indices for strains used in the Mixed ILP and for each sample,
it will output strains infecting the sample and their proportions. 

This script takes 3 arguments: 
1) Path to the directory which contains the information about variants and their proportions at each locus
    for each sample
2) Path to MLST reference file
3) Output path for files 
"""
import pandas as pd
import csv
import os
import re
import itertools
import numpy as np
import cplex
import sys

''' ====================================== Function Definition ======================================================= '''

'''Create a function to read all data files and return a dictionary where keys are the sample names, 
values are dataframes which contain the information about variants and proportions at different loci.
It also return total number of samples and starting sample number(based on last 3 digits)

Input
dataFilePath: File path that contains all your sample folders. dataFilePath should only contain directories of samples and 
              a reference csv
lociOrder: a list that contains order of columns that you want '''
def readData(dataFilePath, lociOrder):
    data = dict()
    sampleFold = list()
    
    for folder in os.listdir(dataFilePath):
        if not folder.endswith(".csv"):
            sampleFold.append(folder)
            
    numSamples = len(sampleFold)
    startingSamp = min([int(i[-3:]) for i in sampleFold])
    
    for folder in sampleFold:
        data["{}".format(folder)] = pd.DataFrame(columns=lociOrder)   #require column to be a specfic order based on lociOrder
        sampleFilePath = dataFilePath + "/{}".format(folder)
        dirs= os.listdir(sampleFilePath)
        csvFiles = [i for i in dirs if i.endswith("proportions.csv")]
        
        temp = re.compile('(.*)_Simulation*')  #grab gene name
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
            
    return data, numSamples, startingSamp

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
            #print column
            temp = (column.keys()[0]).split('_')
            #locus = temp[0]
            #print locus
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
'''
def returnCombinationsAndNumComb(data, numLoci, loci):
    strains = list()
    numOfComb = dict()
    previousNum = 0
    for sample in data:
        #key = sample name, value = dataframe
        sampleDF = data[sample]
        variantsAtAllLoci = list()
        
        #Use for getting matrix representation of combinations
        #maxNumVar = 0
        #only one row in the dataframe
        for column in sampleDF.values[0]:
            variantsAtAllLoci.append(column.keys())
            
            '''
            numVar = len(column.keys())
            if numVar > maxNumVar:
                maxNumVar = numVar
            '''
        
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
    return (strains, numOfComb)

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
            strList = strain[(strain[name]==var) & (strain["Sample"]=="{}".format(sample))]["ST"]
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

def strainSolver(dataPath, refStrains, outputPath, loci):
    ''' ============================================== Data handling ====================================================== '''
    #paramaters
    propFormat = 100    #proportion in percentage or fraction
    #loci = ['clpA', 'clpX', 'nifS']
    numLoci = len(loci)
    
    #read data for samples and reference
    data, numSamples, startingSampleNum = readData(dataPath,loci)
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
    strainAndNumComb = returnCombinationsAndNumComb(data, numLoci, loci)
    strains = strainAndNumComb[0]
    numOfComb = strainAndNumComb[1]
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
    proportionWeightDecVarDF["_merge"] = proportionWeightDecVarDF["_merge"].where(proportionWeightDecVarDF['_merge'] == "left_only", 0)
    proportionWeightDecVarDF["_merge"] = proportionWeightDecVarDF["_merge"].where(proportionWeightDecVarDF['_merge'] == 0, 1)
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
    #model.parameters.mip.pool.intensity.set(2)
    #model.parameters.mip.limits.populate.set(2100000000)
    #model.parameters.mip.pool.absgap.set(0)
    #model.parameters.mip.pool.replace.set(1)
    #model.populate_solution_pool()
   
    objvalue = model.solution.get_objective_value()
    varNames = model.variables.get_names()
    varValues = model.solution.get_values(varNames)
    conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
    conclusion["Decision Variable"] = varNames
    conclusion["Value"] = varValues
    strainInfo = conclusion.merge(strainWeightDecVarDF[strainWeightDecVarDF["Decision Variable"].isin(varNames)])
    strainInfo["New/Existing"] = ["Existing" if w==0 else "New" for w in strainInfo["Weights"].tolist()]
    strainsNeeded = (strainInfo[strainInfo["Value"] > 0][loci + ["ST", "New/Existing"]])
    
    
    #output indices of all strains (New/Existing)
    allStr = strainWeightDecVarDF[["ST", "Weights"] + loci]
    allStr["New/Existing"] = ["Existing" if w==0 else "New" for w in allStr["Weights"].tolist()]
    allStr.drop("Weights", 1, inplace=True)
    allStr.to_csv("{0}/indexedStrains.csv".format(outputPath))
    
    output = proportionWeightDecVarDF[proportionWeightDecVarDF["Sample"] == allSamples[0]].merge(strainsNeeded).drop(["Weights", "Sample"],1)
    output["Proportion"] = model.solution.get_values(output["Decision Variable"].tolist())
    output.drop("Decision Variable", axis=1, inplace=True)
    output = output[loci+["Proportion"]]
    return output
