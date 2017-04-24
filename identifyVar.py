#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:42:23 2017

@author: stanleygan

Integer Linear Program for figuring out variants which cover most of the read

"""
import pandas as pd
import numpy as np
import cplex
import csv
import re
import os
import glob

def returnDataMatrix(path):
    data_matrix = pd.read_csv(path)
    data_matrix.rename(columns={'Unnamed: 0':'Read'}, inplace=True)
    return data_matrix
    
def precision(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(predicted))
        
def recall(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(true))
    
# pred_prop and true_prop: dictionaries
def totalVariationDist(pred_prop, true_prop):
    common_keys = set.intersection(set(pred_prop.keys()), set(true_prop.keys()))
    diff_keys_pred = set.difference(set(pred_prop.keys()), common_keys)
    diff_keys_true = set.difference(set(true_prop.keys()), common_keys)
    common_keys = list(common_keys)
    diff_keys_pred = list(diff_keys_pred)
    diff_keys_true = list(diff_keys_true)
    
    totalVarDist=0
    for key in common_keys:
        totalVarDist += abs(pred_prop[key] - true_prop[key])
        
    for key in diff_keys_pred:
        totalVarDist += pred_prop[key]
        
    for key in diff_keys_true:
        totalVarDist += true_prop[key]
        
    return totalVarDist/2
    
def predictCorrectly(predicted, true):
    return set(predicted) == set(true)  #1 if true, 0 if false
    
#def solver(dataFile, fileName, storeFile):
''' dataMatrix: data frame '''
def solver(dataMatrix):
#    data_matrix = returnDataMatrix("/home/glgan/Documents/Borrelia/data/simData/clpA_7_weighted.csv")
#    data_matrix = returnDataMatrix(dataFile+fileName)
    data_matrix = dataMatrix
    total_read = data_matrix.shape[0]
    total_var = data_matrix.shape[1] - 1  #As first column are reads
    
    xIsVariant = pd.DataFrame(index=[var for var in data_matrix.columns.tolist()[1:]], data=["x{0}".format(i+1) for i in range(total_var)], columns=["Variable"])
    variantName_variable_dict = xIsVariant["Variable"].to_dict()
#    variable_variantName_dict = pd.DataFrame(data=[var for var in data_matrix.columns.tolist()[1:]], index=["x{0}".format(i+1) for i in range(total_var)], columns=["Variant"])["Variant"].to_dict()
    reads = list()
    readAndMM = list()
    
    #Grab unique set of mismatches for each read 
    for row in data_matrix.itertuples():
        read = row[1]
        varMM = list(set(list(row[2:])))
        varMM = [int(mm) for mm in varMM if mm != -1]
        reads.append(read)
        readAndMM.append(varMM)
            
    readAndMMVariable = [[[j,"y{0}_{1}".format(i+1, j)] for j in readAndMM[i]] for i in range(total_read)]
    yIsRead = pd.DataFrame(data=reads, columns=["Read"])
    yIsRead["Variable"] = readAndMMVariable
    
    #Dictionary for easier retrievable of y_variable name when formulating ILP
    #key is a tuple of read name and number of mismatches, value = read variable name
    readMM_varName_dict = dict()
    #key is y variable name, value is the matching read
    varName_read_dict = dict()
    for r in range(total_read):
        for mm in readAndMM[r]:
            readMM_varName_dict[(reads[r], mm)] = "y{0}_{1}".format(r+1, mm)
            varName_read_dict["y{0}_{1}".format(r+1, mm)] = reads[r]
            
           
    '''=============================================== Form ILP ====================================================== '''
    model = cplex.Cplex()
    model.objective.set_sense(model.objective.sense.minimize)
    #add variables related to variants, x_j represents variant j
    model.variables.add(obj=[1]*total_var, names=xIsVariant["Variable"].tolist(), types=[model.variables.type.binary]*total_var)
    #model.variables.add(names=xIsVariant["Variable"].tolist(), types=[model.variables.type.binary]*total_var)
    
    #add variables related to reads, y_ik means read i with k mismatches
    y_weights = [mm[0] for mmName in readAndMMVariable for mm in mmName]
    y_variables = [mm[1] for mmName in readAndMMVariable for mm in mmName]
    model.variables.add(obj=y_weights, names=y_variables, types=[model.variables.type.binary]*len(y_variables))
    
    #Constraint: Cover most of the read
    keep_alpha = 0.99
    model.linear_constraints.add(lin_expr=[[y_variables, [1]*len(y_variables)]], senses=["G"], rhs=[keep_alpha*total_read], names=["c{0}".format(1+model.linear_constraints.get_num())])
    
    #Constraint: For each read i, only allow it to be covered with a unique number of mismatch i.e. 
    #only one of y_i0, y_i1, y_i2, y_i3 is =1
    uniqueMMConstr = [[[mmName[1] for mmName in oneRead], [1]*len(mmName[1])] for oneRead in readAndMMVariable ]
    model.linear_constraints.add(lin_expr=uniqueMMConstr, senses=["L"]*total_read, rhs=[1]*total_read, names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(total_read)])
    
    #Limit maximum number of variants
    #model.linear_constraints.add(lin_expr=[[xIsVariant["Variable"].tolist(), [1]*total_var]], senses=["L"], rhs=[10], names=["c{0}".format(1+model.linear_constraints.get_num())])
    
    #This is for easier formulation of the constraints related to: If y_ik=1, then it must be 
    #covered by some variant with k mm
    yCoverConstr = list()
    readVar_variant_dict = dict()   #key is read variable name, value is the list of variants covering that read(with the particular mm)
    for row in data_matrix.itertuples():
        read = row[1]
        uniqueMm = list(set(row[2:]))
        uniqueMm = [int(mm) for mm in uniqueMm if mm!=-1]
        
        for mm in uniqueMm:
            temp = [i for i in range(len(list(row[1:]))) if list(row[1:])[i] == mm]
            
            xVariable = [variantName_variable_dict[variant] for variant in data_matrix.columns[temp].tolist()]
            yVariable = [ readMM_varName_dict[read, mm] ]
            readVar_variant_dict[readMM_varName_dict[read, mm]] = xVariable
            yCoverConstr.append([ xVariable + yVariable, [1]*len(xVariable) + [-1]*len(yVariable) ])
            
    #Constraint: If y_ik is chosen, it must be covered by some variant j with k mm
    model.linear_constraints.add(lin_expr=yCoverConstr, senses=["G"]*len(y_variables), rhs=[0]*len(y_variables), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(y_variables))])
    
    #model.write("a.lp")
    model.set_results_stream(None)
    model.solve()
       
    '''=========================================== Presenting Results ========================================================'''
    objvalue = model.solution.get_objective_value()
    varNames = model.variables.get_names()
    varValues = model.solution.get_values(varNames)
    conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
    conclusion["Decision Variable"] = varNames
    conclusion["Value"] = varValues

    present = conclusion[conclusion["Value"]==1]
    variantsPresent = xIsVariant[xIsVariant["Variable"] .isin(present["Decision Variable"].tolist())]
    varPresentList = variantsPresent.index.tolist()
#    xPresentList = variantsPresent["Variable"].tolist()

    yCovered = present[present["Decision Variable"].str.contains(u"y.*")]["Decision Variable"].tolist()
    readsCovered = [varName_read_dict[y] for y in yCovered]
    return objvalue, varPresentList, readsCovered
    
'''    
    #Output predicted variants
    name = re.sub("_weighted.csv", "" ,fileName)
    with open(storeFile+name+"_variants.csv", "w") as f:
        writer = csv.writer(f)
        for var in variantsPresent.index.tolist():
            writer.writerow([var])

    #Find assignment of reads to variants and the proportions of the variants
    predictedProp = dict.fromkeys(variantsPresent["Variable"].tolist(), 0.0)
    readsCovered = present[present["Decision Variable"].str.contains(u"y.*")]["Decision Variable"].tolist()
    
    for readVar in readsCovered:
        variantCovering = readVar_variant_dict[readVar]
        variantCoveringAndPredicted = list(set.intersection(set(variantCovering), set(xPresentList)))
        
        for variant in variantCoveringAndPredicted:
            increment = 1/len(variantCoveringAndPredicted)
            predictedProp[variant] = predictedProp[variant] + increment
                         
    #Normalize the proportions
    normalize_term = 1.0/sum(predictedProp.values())
    
    for key in predictedProp:
        predictedProp[key] = normalize_term * predictedProp[key]
        
    #Match x variables back to variant name
    variantName_proportions = dict()
    for key in predictedProp:
        variantName_proportions[variable_variantName_dict[key]] = predictedProp[key]*100
        
    #output the predicted proportions
    with open(storeFile+name+"_predicted_proportions.csv", "w") as f:
        writer = csv.writer(f)
        for key, val in variantName_proportions.items():
            writer.writerow([key,val])
'''
 
''' ================================= Main ============================================= '''
if __name__ == "__main__":
    #identifyVar()
    dataFile = "/home/glgan/Documents/Borrelia/data/simData/"
    storeFile = "/home/glgan/Documents/Borrelia/data/predictedVar/"
    #dataFile = "/home/glgan/Documents/borrelia/data/preproc/"
    #fileName = "toy.csv"
    extension = "csv"
    os.chdir(dataFile)
    csvFiles = [f for f in glob.glob("*.{}".format(extension))]
    
    epoch=0
    for fileName in csvFiles:
        print("Epoch:{} ".format(epoch))
        solver(dataFile, fileName, storeFile)
        epoch +=1
