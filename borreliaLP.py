#!/usr/bin/env python2
# -*- coding: utf-8 -*-
''' Remark: Have to check correct proportions after reading file, convert proportions to 3 dec. 
Transfer zip to itertools'''
"""
Created on Mon Jan 16 11:33:50 2017

@author: stanleygan
"""
import pandas as pd
import csv
import os
import re
import itertools
import numpy as np
import cplex

''' ====================================== Function Definition ======================================================= '''

'''Create a function to read all data files and return a dictionary where keys are the sample names, 
values are dataframes which contain the information about variants and proportions at different loci

Input
dataFilePath: File path that contains all your sample folders 
numSamples: number of samples
startingSampleNum: the unique index of the first sample
lociOrder: a list that contains order of columns that you want '''
def readData(dataFilePath, numSamples, startingSampleNum, lociOrder):
    data = dict()
    
    for i in range(startingSampleNum, startingSampleNum+numSamples):
        data['sample%d' %i] = pd.DataFrame(columns=lociOrder)   #require column to be a specfic order based on lociOrder
        
        sampleFilePath = dataFilePath + ('/sample%d/' %i)
        dirs= os.listdir(sampleFilePath)
        for files in dirs:
            temp = re.compile('_(.*).csv')  #grab gene name
            gene = temp.findall(files)
            
            reader = csv.reader(open(sampleFilePath+files, 'r'))
            
            #store variants and respective proportions in a dictionary
            d = dict()  
            for variant, proportion in reader:
                d[variant] = proportion
            
            #wrap dictionary with a list for storage in dataframe
            alist =[d]
            (data['sample%d' %i])[gene[0]] = alist
            
    return data

def roundAndCheckProp(data):
    for sample in data:
        sampleDF = data[sample]
        
        for column in sampleDF.values[0]:
            prop = column.values()
            prop = [float(p) for p in prop]
            
   
#'''Return the unique combinations of variants at all loci across all samples
#data: Data file preloaded previously'''
#def uniqueCombinations(data):
#    uniqueStrains = list()
#    for sample in data:
#        #key = sample name, value = dataframe
#        sampleDF = data[sample]
#        variantsAtAllLoci = list()
#        #only one row in the dataframe
#        for column in sampleDF.values[0]:
#            variantsAtAllLoci.append(column.keys())
#        
#        combination = itertools.product(*variantsAtAllLoci)
#        for strain in combination:
#            uniqueStrains.append(strain)
#            
#    uniqueStrains = list(set(uniqueStrains))
#    return uniqueStrains

'''Returns a dataframe of combinations of the variants at all loci, the matrix representation of the combination
and the sample that contains it. Also, it returns the total number of combinations that each sample has

Input
data: dictionary, information of all samples preloaded previously
numLoci: number of loci
'''
def returnCombinationsAndNumComb(data, numLoci):
    strains = list()
    numOfComb = dict()
    previousNum = 0
    for sample in data:
        #key = sample name, value = dataframe
        sampleDF = data[sample]
        variantsAtAllLoci = list()
        maxNumVar = 0
        #only one row in the dataframe
        for column in sampleDF.values[0]:
            variantsAtAllLoci.append(column.keys())
            
            numVar = len(column.keys())
            if numVar > maxNumVar:
                maxNumVar = numVar
        
        combination = itertools.product(*variantsAtAllLoci)
        combinationIndices = [list(comb) for comb in itertools.product(*[range(len(var)) for var in variantsAtAllLoci])]
        for strain,strainIndex in itertools.izip(combination,combinationIndices):
            temp = list(strain)
            temp.append(sample) #add sample name
            matrixRep = np.zeros(shape=(maxNumVar, numLoci))
            matrixRep[strainIndex, np.arange(numLoci)] = 1
            temp.append(matrixRep)
            strains.append(temp)
            
        #Get the total number of combinations for each sample
        numOfComb[sample] = len(strains) - previousNum
        previousNum = numOfComb[sample]
        
    strains = pd.DataFrame(strains, columns=(loci+['Sample', 'Matrix Representation']))
    return (strains, numOfComb)

'''Returns a dictionary where keys=sample name and values=numpy matrix representation of the proportions for the sample
Input
data: dictionary, information of all samples preloaded previously
numLoci: number of loci
'''
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
        d = data[sample]
        for column in d:
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
'''                                                                           
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
Return a data frame with locus as columns, and each column contains a dictionary(wrapped by a list) which shows the indices of unique strains(ST) that each key=(variants, sample) maps to
Input:
    strain: a dataframe with variants combination
    loci : a list containing locus name
    start: starting sample number
    numSamp: number of samples
'''                                                                           
def mapVarAndSampleToStrain(strain, loci, start, numSamp):
    varToStr = pd.DataFrame(columns=loci)
    for name in loci:
        varDF = strain.drop_duplicates(name)
        varDF = (varDF[name]).reset_index(drop=True)
        varDict = dict()
        for sample, var in list(itertools.product(range(start, start+numSamp), varDF.tolist())):
            strList = strain[(strain[name]==var) & (strain["Sample"]=="sample{0}".format(sample))]["ST"]
            if len(strList) != 0:
                varDict[(var, sample)] = strList
               
        varToStr[name] = [varDict]
        
    return varToStr




''' ============================================== Data handling ====================================================== '''
#paramaters
numLoci = 3
startingSampleNum = 1
numSamples = 2
loci = ['clpA', 'clpX', 'nifS']

#read data for samples and reference
data = readData("/home/stanleygan/Documents/Borrelia/data/simpleEx",numSamples,startingSampleNum, loci)
reference = pd.read_csv('~/Documents/Borrelia/data/simpleEx/reference.csv',usecols=range(1,numLoci+1))
lociNames = list(reference.columns.values)
numReference = reference.shape[0]

roundAndCheckProp(data)

#As reference only contains numbers as entries, add gene name to the variants for better identification
for name in lociNames:
    reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
    
#Get proportions of variants at different locus for each sample
varAndProp = returnVarAndProportions(data)

#Get the combinations at all loci across all samples
strainAndNumComb = returnCombinationsAndNumComb(data, numLoci)
strains = strainAndNumComb[0]
numOfComb = strainAndNumComb[1]
uniqueStrains = strains.drop_duplicates(loci)
uniqueStrains = (uniqueStrains[loci]).reset_index(drop=True)
uniqueStrains["ST"] = uniqueStrains.index.values + 1    #assign indices for strains or each unique combinations
strains = strains.merge(uniqueStrains, indicator=True, how="left")    #assign the index to the combinations(as strain data frame contains duplicated rows)
strains = strains.drop("_merge",1)

#For each variants, get a mapping of which strains it maps to
varSampToST = mapVarAndSampleToStrain(strains, loci, startingSampleNum, numSamples)
#for a in varSampToST:
#    dictionary = varSampToST[a][0]
#mapVtoS = mapVarToStrain(uniqueStrains, loci)
#clpa = mapVtoS["clpA"][0]

#weights and decision variables for proportion of strains. weight=0 if the strain is in reference, otherwise =1. Notice there will be duplications of strain types
#here because for proportions, we consider sample by sample rather than unique strain types
proportionWeightDecVarDF = strains.merge(reference, indicator=True, how="left")
proportionWeightDecVarDF["_merge"] = proportionWeightDecVarDF["_merge"].where(proportionWeightDecVarDF['_merge'] == "left_only", 0)
proportionWeightDecVarDF["_merge"] = proportionWeightDecVarDF["_merge"].where(proportionWeightDecVarDF['_merge'] == 0, 1)
proportionWeightDecVarDF = proportionWeightDecVarDF.rename(columns = {"_merge":"Weights"})

#Add proportion decision variable names
proportionWeightDecVarDF["Decision Variable"] = np.nan
for i in range(startingSampleNum, startingSampleNum + numSamples):
    propNameTemp = ["pi_%s_%d" %t for t in itertools.izip((proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == "sample%d" %i])['Sample'], range(1,1+(proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == "sample%d" %i]).shape[0]))]
    propNameTemp = [ele.replace("pi_sample", "pi_s") for ele in propNameTemp]   #shorter name as CPLEX can't hold name with >16 char
    proportionWeightDecVarDF.loc[proportionWeightDecVarDF['Sample'] == "sample%d" %i, 'Decision Variable'] = propNameTemp
    
#weights and decision variables for unique strain types, weight=0 if strain is in reference, otherwise=1. no duplications
strainWeightDecVarDF = proportionWeightDecVarDF.drop_duplicates(loci)
retainCol = loci + ['Weights', 'ST']
strainWeightDecVarDF = strainWeightDecVarDF[retainCol].reset_index(drop=True)
strainWeightDecVarDF["Decision Variable"] = ["a%d" %(i+1) for i in range(len(strainWeightDecVarDF['Weights'].values.tolist()))]

'''==================================== Forming LP here ================================================'''

#Form a CPLEX model
model = cplex.Cplex()
#minimize problem
model.objective.set_sense(model.objective.sense.minimize)
#add the decision variables for unqiue strain types
model.variables.add(obj=strainWeightDecVarDF['Weights'].values.tolist(), names=strainWeightDecVarDF['Decision Variable'], types = [model.variables.type.binary]* len(strainWeightDecVarDF['Weights'].values.tolist()))
#add proportions decision variables
model.variables.add(obj=proportionWeightDecVarDF['Weights'].values.tolist(),ub=[1]*proportionWeightDecVarDF['Weights'].shape[0], names=proportionWeightDecVarDF["Decision Variable"], types=[model.variables.type.continuous] * len(proportionWeightDecVarDF['Weights'].values.tolist()))

#add linear constraints such that for each sample, the sum of the proportions of its variants combination = 1
propVarSumTo1 = list()

for i in range(numSamples):
    temp=list()
    for num in range(1, 1+ numOfComb["sample{0}".format(i + startingSampleNum)]):
        temp.append("pi_s{0}_{1}".format(i + startingSampleNum, num))
        
    propVarSumTo1.append([temp, [1]* len(temp)])
    
model.linear_constraints.add(lin_expr=propVarSumTo1, rhs=[1]*len(propVarSumTo1), senses=["E"]*len(propVarSumTo1), names=["c{0}".format(i+1) for i in range(len(propVarSumTo1))])

#add linear constraints such that for each sample, sum of pi_ik \dot V_ik (proportion \dot matrix representation) across all combinations = Proportion matrix
piDotComb = list()
propConstrRHS = list()
for name in varSampToST:
    temp=list()
    varSampToSTDict = varSampToST[name][0]
    for (var, sample) in varSampToSTDict:
        strainTypes = varSampToSTDict[(var, sample)]
        propDecVar = proportionWeightDecVarDF[(proportionWeightDecVarDF["ST"].isin(strainTypes)) & (proportionWeightDecVarDF["Sample"] == "sample{0}".format(sample))]["Decision Variable"]
        propConstrRHS.append(  float( ( (data["sample{0}".format(sample)])[name][0] )[var] )  )
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
    indicMinusProp.append([[i, pi],[1, -1]])    

model.linear_constraints.add(lin_expr=indicMinusProp, rhs=[0]*len(indicMinusProp), senses=["G"]*len(indicMinusProp), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(indicMinusProp))] )

#add error variables and linear constraints related to error terms
#create error variable names
errorName = ["d_"]*(varAndProp.shape[0]) + varAndProp["Sample"] + ["_"]*(varAndProp.shape[0]) + varAndProp["Variant"]
errorName = [ele.replace("d_sample", "d_s") for ele in errorName]
varAndProp["Error"] = errorName
          
#add error variable
model.variables.add(obj=[1]*varAndProp.shape[0], lb=(-1*varAndProp["Proportion"]).tolist(), ub=(1-varAndProp["Proportion"]).tolist(), names=varAndProp["Error"].tolist(), types=[model.variables.type.continuous]*varAndProp.shape[0])

#add the constraints whereby for each sample, at each locus, the sum of the error of all variants=0
errorSumTo0 = list()
for samp, loc in list(set(itertools.izip(varAndProp["Sample"].tolist(), varAndProp["Locus"].tolist()))):
    temp = (varAndProp[(varAndProp["Sample"] == samp) & (varAndProp["Locus"] == loc)])["Error"].tolist()
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
    
    err = row["Error"]
    pi = proportionWeightDecVarDF[( proportionWeightDecVarDF["Sample"]  == samp ) & (proportionWeightDecVarDF[loc] == var )]["Decision Variable"].tolist()
    prop = row["Proportion"]
    
    errLessSumMinProp.append( [[err] + pi, [i for i in itertools.chain([1],[-1]*len(pi))]] )
    errLessSumMinPropRHS.append(-1*prop)
    
    errLessPropMinSum.append( [[err] + pi, [i for i in itertools.chain([1],[1]*len(pi))]] )
    errLessPropMinSumRHS.append(prop)

model.linear_constraints.add(lin_expr=errLessSumMinProp, rhs=errLessSumMinPropRHS, senses=["L"]*len(errLessSumMinProp), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errLessSumMinProp))])  
model.linear_constraints.add(lin_expr=errLessPropMinSum, rhs=errLessPropMinSumRHS, senses=["L"]*len(errLessPropMinSum), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(errLessPropMinSum))])

#Final part
#model.write("a.lp")
#model.solve()
#objvalue = model.solution.get_objective_value()
#varNames = model.variables.get_names()
#varValues = model.solution.get_values(varNames)
#conclusion = pd.DataFrame(columns=["Name", "Value"])
#conclusion["Name"] = varNames
#conclusion["Value"] = varValues

#varNames = [name for name in model.variables.get_names() if any(char in name for char in 'pi_')]
#varValues = model.solution.get_values(varNames)
#conclusion = pd.DataFrame(columns=["Name", "Value", "Existing(0) or New(1)"])
#conclusion["Name"] = varNames
#conclusion["Value"] = varValues
#conclusion[["Existing(0) or New(1)", "ST"] + loci] = proportionWeightDecVarDF[proportionWeightDecVarDF["Decision Variables"].isin(varNames)][["Weights", "ST"] + loci]
#newStrainIntro = conclusion[(conclusion["Existing(0) or New(1)"] == 1) & (conclusion["Value"] > 0)]