#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:50:06 2017

@author: stanleygan

Integer Linear Program for figuring out variants which cover most of the read

"""
import pandas as pd
import numpy as np
import cplex
import itertools

class preprocessing:
    def __init__(self, sample, locus, path):
        self.data_matrix = pd.read_csv(path)
        self.data_matrix.rename(columns={'Unnamed: 0':'Read'}, inplace=True)
        self.sample = sample
        self.locus = locus
        self.total_read = self.data_matrix.shape[0]
        self.total_var = self.data_matrix.shape[1] - 1  #As first column are reads
        #key = variant, value=list of reads
        self.readsForEachVar = dict()
        #key = read, value=list of variants covering it
        self.varsForEachRead = dict()
        
        #some decision variables
        self.xIsVariant, self.yIsRead = self._setVariables()
        
        #handle some data
        self._getReadsForEachVar()
        self._getVarsForEachRead()
        
    def _setVariables(self):
        xIsVariant = pd.DataFrame(data=[var for var in self.data_matrix.columns.tolist()[1:]], index=["x{0}".format(i+1) for i in range(self.total_var)], columns=["Variant"])
        yIsRead = pd.DataFrame(data=[read for read in self.data_matrix["Read"].tolist()], index=["y{0}".format(i+1) for i in range(self.total_read)], columns=["Read"])
        
        return xIsVariant, yIsRead
    
    def _getReadsForEachVar(self):
        for var in self.data_matrix.columns.tolist()[1:]:
            reads = self.data_matrix[self.data_matrix[var] == "1.0"]["Read"].tolist()
            variables = self.yIsRead[self.yIsRead["Read"].isin(reads)].index.tolist()
            variableRespToVar = "".join(self.xIsVariant[self.xIsVariant["Variant"] == var].index.tolist())
            self.readsForEachVar[variableRespToVar] = variables
                                 
    def _getVarsForEachRead(self):
        for read in self.data_matrix["Read"].tolist():
            temp = self.data_matrix[self.data_matrix["Read"] == read ] == "1.0"
            var = self.data_matrix.columns[np.where(temp.iloc[0,:] == True)]
            variables = self.xIsVariant[self.xIsVariant["Variant"].isin(var)].index.tolist()
            variableRespToRead = "".join(self.yIsRead[self.yIsRead["Read"] == read].index.tolist())
            self.varsForEachRead[variableRespToRead] = variables
            
    def getter_readsForEachVar(self):
        return self.readsForEachVar
    
    def getter_varsForEachRead(self):
        return self.varsForEachRead
    
    def getter_dataMatrix(self):
        return self.data_matrix
    
    def getter_totalRead(self):
        return self.total_read
    
    def getter_totalVar(self):
        return self.total_var
    
    def getter_decisionVariables(self):
        return self.xIsVariant, self.yIsRead
        
if __name__ == "__main__":
    p = preprocessing("Sample1", "clpA", "~/Downloads/clpA_7.csv")
    a = p.getter_dataMatrix()
    readsForEachVar = p.getter_readsForEachVar()
    varsForEachRead = p.getter_varsForEachRead()
    totalVar = p.getter_totalVar()
    totalRead = p.getter_totalRead()
    xIsVariant, yIsRead = p.getter_decisionVariables()
    keep_alpha = 0.99
    
    ''' ============================ Build ILP model here ============================ '''
    #Idea is to minimize the number of variants while covering most reads
    
    #create a mapping from variables to name of variants/reads
    model = cplex.Cplex()
    model.objective.set_sense(model.objective.sense.minimize)
    #add variables related to variants, x_i represents variant i
    model.variables.add(obj=[1]*totalVar, names=xIsVariant.index.tolist(), types=[model.variables.type.binary]*totalVar)
    #add variables related to reads, y_j means read j
    model.variables.add(names=yIsRead.index.tolist(), types=[model.variables.type.binary]*totalRead)
    
    #Constraint:Most reads should be covered
    model.linear_constraints.add(lin_expr=[[yIsRead.index.tolist(), [1]*totalRead]], senses=["G"], rhs=[keep_alpha*totalRead], names=["c{0}".format(1+model.linear_constraints.get_num())])
    
    #Constraint: If a read is chosen, then it must be covered by some variant
    readCoverConstr = [[varsForEachRead[key] + [key], [1]*len(varsForEachRead[key]) + [-1]] for key in varsForEachRead.keys()]
    model.linear_constraints.add(lin_expr=readCoverConstr, senses=["G"]*len(readCoverConstr), rhs=[0]*len(readCoverConstr), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(readCoverConstr))])
    #model.write("hello.lp")
    
    #solve
    model.solve()
    objvalue = model.solution.get_objective_value()
    varNames = model.variables.get_names()
    varValues = model.solution.get_values(varNames)
    conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
    conclusion["Decision Variable"] = varNames
    conclusion["Value"] = varValues
              
    present = conclusion[conclusion["Value"]==1]