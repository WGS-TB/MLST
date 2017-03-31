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

def returnDataMatrix(path):
    data_matrix = pd.read_csv(path)
    data_matrix.rename(columns={'Unnamed: 0':'Read'}, inplace=True)
    return data_matrix

data_matrix = returnDataMatrix("/home/stanleygan/Documents/Borrelia/data/preproc/toy.csv")
total_read = data_matrix.shape[0]
total_var = data_matrix.shape[1] - 1  #As first column are reads

xIsVariant = pd.DataFrame(index=[var for var in data_matrix.columns.tolist()[1:]], data=["x{0}".format(i+1) for i in range(total_var)], columns=["Variable"])
xIsVariantDict = xIsVariant["Variable"].to_dict()
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
yIsReadDict = dict()
for r in range(total_read):
    for mm in readAndMM[r]:
        yIsReadDict[(reads[r], mm)] = "y{0}_{1}".format(r+1, mm)
        
       
'''====================================== Form ILP ============================================= '''
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
model.linear_constraints.add(lin_expr=[[xIsVariant["Variable"].tolist(), [1]*total_var]], senses=["L"], rhs=[10], names=["c{0}".format(1+model.linear_constraints.get_num())])

#This is for easier formulation of the constraints related to: If y_ik=1, then it must be 
#covered by some variant with k mm
yCoverConstr = list()
for row in data_matrix.itertuples():
    read = row[1]
    uniqueMm = list(set(row[2:]))
    uniqueMm = [int(mm) for mm in uniqueMm if mm!=-1]
    
    for mm in uniqueMm:
        temp = [i for i in range(len(list(row[1:]))) if list(row[1:])[i] == mm]
        
        xVariable = [xIsVariantDict[variant] for variant in data_matrix.columns[temp].tolist()]
        yVariable = [ yIsReadDict[read, mm] ]
        yCoverConstr.append([ xVariable + yVariable, [1]*len(xVariable) + [-1]*len(yVariable) ])
        
#Constraint: If y_ik is chosen, it must be covered by some variant j with k mm
model.linear_constraints.add(lin_expr=yCoverConstr, senses=["G"]*len(y_variables), rhs=[0]*len(y_variables), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(y_variables))])

#model.write("a.lp")
model.solve()

objvalue = model.solution.get_objective_value()
varNames = model.variables.get_names()
varValues = model.solution.get_values(varNames)
conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
conclusion["Decision Variable"] = varNames
conclusion["Value"] = varValues
          
present = conclusion[conclusion["Value"]==1]
        
#if __name__ == "__main__":
    #identifyVar()