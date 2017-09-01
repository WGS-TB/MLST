import os
import getVariantAndProp as gvp
import getStrAndProp as gsp
import numpy as np
import itertools
import pandas as pd
import csv
import time
import argparse
import cplex

start_time = time.time()
ap = argparse.ArgumentParser()
ap.add_argument("-s", "--sample", required = False, default="all")
args = vars(ap.parse_args())

#currentpath = /pipeline/
currentPath = os.getcwd()
data_path = currentPath +"/data/"
lociDb_path = currentPath + "/loci_db/"
ref_strains = currentPath + "/strain_ref.txt"
loci = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]
reference = pd.read_csv(currentPath+"/strain_ref.txt",sep="\t",usecols=range(1,len(loci)+1))
for name in loci:
    reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)

if not os.path.exists("strainsAndProp_errThres"):
    os.mkdir("strainsAndProp_errThres")

if args["sample"] != "all":
    gsp.strainSolver(currentPath+"/variantsAndProp/{}".format(args["sample"]), currentPath+"/strain_ref.txt", currentPath+"/strainsAndProp_errThres", loci, "all", args["sample"])
else:
    gsp.strainSolver(currentPath+"/variantsAndProp", currentPath+"/strain_ref.txt", currentPath+"/strainsAndProp_errThres", loci, "all", args["sample"])
#if args["sample"] != "all":
#    dataPath = currentPath+"/variantsAndProp/{}".format(args["sample"])
#else:
#    dataPath = currentPath+"/variantsAndProp"
#
#print("\n~~~~~~~~~~~~~~~~~~~~~~ Solving ILP for strain prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
#solution_dict, ilp_objective_dict, data, strains = gsp.minNewStrain(dataPath, currentPath+"/strain_ref.txt", loci, args["sample"])
##    print(solution_dict)
#print("Number of solutions from ILP: {}".format(len(solution_dict)))
#print("\n~~~~~~~~~~~~~~~~~~~~~~ Solving LP for proportion prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
#feasible_sol = list()
#lp_objective_dict = dict()
#sample_strainProp_dict = dict()
#infeasibility = 0
#for i in solution_dict.keys():
#    print("============= ILP Solution {0} ===============".format(i))
#    feasible=False
#    try:
#        objvalue, sampleAndStrainProp = gsp.minNewStrainProp(solution_dict[i], data, strains, ref_strains, loci)
#        feasible=True
#    except cplex.exceptions.errors.CplexSolverError as e:
#        infeasibility +=1
#        print("Infeasible!")
#        
#    if feasible == True:
#        print("Objective Value (Pure ILP, LP, Sum): ({0}, {1}, {2})".format(ilp_objective_dict[i], objvalue, objvalue+ilp_objective_dict[i]))
#        feasible_sol.append(i)
#        lp_objective_dict[i] = objvalue
#        #sampleAndStrainProp is a dict with samples as keys and dataframes as values
#        sample_strainProp_dict[i] = sampleAndStrainProp
#    
#if infeasibility == len(solution_dict):
#    print("No feasible solutions")
#else:
#    total_obj_dict = dict()
#    min_obj = np.inf
#    for j in feasible_sol:
#        total_obj_dict[j] = ilp_objective_dict[j] + lp_objective_dict[j]
#        if total_obj_dict[j] < min_obj:
#            min_obj = total_obj_dict[j]
#    
#    min_obj_list = [i for i in total_obj_dict.keys() if total_obj_dict[i] == min_obj]
#    print("\nMinimum objective value: {}".format(min_obj))
#    print("\nSolution which minimizes both pure ILP and LP: {}\n".format(min_obj_list))
#    if len(min_obj_list) > 1:
#        print("@@@@@@@@@@@@@@@ More than 1 optimal solution for the problem @@@@@@@@@@@@@@@@@@")
#        
#    for samp in sample_strainProp_dict[min_obj_list[0]].keys():
#        output = sample_strainProp_dict[min_obj_list[0]][samp]
#        print output
#        output.to_csv("{0}/{1}_strainsAndProportions.csv".format(currentPath+"/strainsAndProp_sep_ALL", samp))
    
print("")
print("Script done.")
print("Created folder strainsAndProp which contains strains and their proportions for each sample")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
