#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 11:34:31 2017

@author: glgan

"""

import pandas as pd
import numpy as np
import cplex


''' dataMatrix: data frame
    Return: Objective value, variants predicted, reads covered by these variants, all optimal solutions, objective values of optimal solutions
'''
def solver(dataMatrix):
#    data_matrix = returnDataMatrix("/home/glgan/Documents/Borrelia/data/simData/clpA_7_weighted.csv")
#    data_matrix = returnDataMatrix(dataFile+fileName)
    data_matrix = dataMatrix
    total_read = data_matrix.shape[0]
    total_var = data_matrix.shape[1]  #As first column are reads

    xIsVariant = pd.DataFrame(index=[var for var in data_matrix.columns.tolist()], data=["x{0}".format(i+1) for i in range(total_var)], columns=["Variable"])
    variantName_variable_dict = xIsVariant["Variable"].to_dict()
#    variable_variantName_dict = pd.DataFrame(data=[var for var in data_matrix.columns.tolist()[1:]], index=["x{0}".format(i+1) for i in range(total_var)], columns=["Variant"])["Variant"].to_dict()
    reads = data_matrix.index.tolist()
    readAndMM = list()

    #Grab unique set of mismatches for each read
    for row in data_matrix.itertuples():
        varMM = list(set(list(row[1:])))
        varMM = [int(mm) for mm in varMM if mm != -1]
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
    keep_alpha = 1.0
    model.linear_constraints.add(lin_expr=[[y_variables, [1]*len(y_variables)]], senses=["G"], rhs=[keep_alpha*total_read], names=["c{0}".format(1+model.linear_constraints.get_num())])

    #Constraint: For each read i, only allow it to be covered with a unique number of mismatch i.e.
    #only one of y_i0, y_i1, y_i2, y_i3 is =1
    uniqueMMConstr = list()
    for oneRead in readAndMMVariable:
        temp = list()
        yVariable = list()
        for i in range(len(oneRead)):
            yVariable.append( (oneRead[i])[1] )
        temp.append(yVariable)
        temp.append([1]*len(yVariable))
        uniqueMMConstr.append(temp)

    model.linear_constraints.add(lin_expr=uniqueMMConstr, senses=["L"]*total_read, rhs=[1]*total_read, names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(total_read)])

    #Limit maximum number of variants
    #model.linear_constraints.add(lin_expr=[[xIsVariant["Variable"].tolist(), [1]*total_var]], senses=["L"], rhs=[10], names=["c{0}".format(1+model.linear_constraints.get_num())])

    #This is for easier formulation of the constraints related to: If y_ik=1, then it must be
    #covered by some variant with k mm
    yCoverConstr = list()
    readVar_variant_dict = dict()   #key is read variable name, value is the list of variants covering that read(with the particular mm)
    for row in data_matrix.itertuples():
        read = row[0]
        uniqueMm = list(set(row[1:]))
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
#    model.solve()

    #options for searching more optimal solutions
    #model.parameters.mip.pool.capacity.set(10)
    model.parameters.mip.pool.intensity.set(4)
    #model.parameters.mip.limits.populate.set(2100000000)
    model.parameters.mip.pool.absgap.set(0)
    model.parameters.mip.pool.replace.set(1)
    model.populate_solution_pool()

    '''=========================================== Presenting Results ========================================================'''
    objvalue = model.solution.get_objective_value()
    varNames = model.variables.get_names()
    varValues = model.solution.get_values(varNames)
    conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
    conclusion["Decision Variable"] = varNames
    conclusion["Value"] = varValues

    #Python numerical instability, sometimes comparison to 1.0 gives false when the value is indeed 1. As all decision variables have max value =1, hence setting > 0.5 is the same as ==1
    present = conclusion[conclusion["Value"]> 0.5]
    variantsPresent = xIsVariant[xIsVariant["Variable"] .isin(present["Decision Variable"].tolist())]
    varPresentList = variantsPresent.index.tolist()
#    xPresentList = variantsPresent["Variable"].tolist()

    yCovered = present[present["Decision Variable"].str.contains(u"y.*")]["Decision Variable"].tolist()
    readsCovered = [varName_read_dict[y] for y in yCovered]

    allSol = list()
    allObjValue = list()
    for i in range(model.solution.pool.get_num()):
        objvalue = model.solution.pool.get_objective_value(i)
        varNames = model.variables.get_names()
        varValues = model.solution.pool.get_values(i,varNames)
        conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
        conclusion["Decision Variable"] = varNames
        conclusion["Value"] = varValues

        present = conclusion[conclusion["Value"]>0.5]
        variantsPresent = xIsVariant[xIsVariant["Variable"] .isin(present["Decision Variable"].tolist())]
        varPresentList = variantsPresent.index.tolist()
        allSol.append(varPresentList)
        allObjValue.append(objvalue)

    return objvalue, varPresentList, readsCovered, allSol, allObjValue
