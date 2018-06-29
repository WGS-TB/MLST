import os
import pipeline_functions as pf
import numpy as np
import pandas as pd
import time
import argparse
import cplex

start_time = time.time()
ap = argparse.ArgumentParser()
#ap.add_argument("-s", "--sample", required = False, default="all", help="Sample name. Default is all samples")
ap.add_argument("-o", "--output", required=False, help="Name of output folder. Default is strainsAndProp", default="strainsAndProp")
#ap.add_argument("-go", "--globalOption", required=True, help="Version of optimization program to use. 'mixed': mixed ILP, 'separated': pure ILP + LP")
#Only for mixed version
ap.add_argument("-oc", "--objectiveComponent", required=False, default="all", help="Objective components. Default: 'all'. 'noProp': Does not include proportion component, 'noPropAndErr': Does not include proportion and error in objective function")
ap.add_argument("-timelim", "--timeLimit", required=False, help="Time limit in integer(sec) for cplex solver for mixed ILP. Default:600sec", default=600)
ap.add_argument("-g", "--gap", required=False, help="Relative gap tolerance(in percent) for cplex solver for mixed ILP. Default: 5", default=5)
ap.add_argument("-pathToDistMat", "--ptdm", required=True, help="Folder name which contains editDistance matrices")
#ap.add_argument("-r", "--ref", required=False, help="Reference strains file name", default="strain_ref.txt")

args = vars(ap.parse_args())
globalSamp = "all"

#currentpath = /pipeline/
currentPath = os.getcwd()
ref_strains = os.path.join(currentPath, "strain_ref.txt")
#ref_strains = currentPath + "/" + args["ref"]
loci = ["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"]
reference = pd.read_csv(ref_strains,sep="\t",usecols=range(1,len(loci)+1))
for name in loci:
    reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)

if not os.path.exists(os.path.abspath(args["output"])):
    os.mkdir(os.path.abspath(args["output"]))

pf.strainSolver(os.path.join(currentPath,"variantsAndProp"), ref_strains, os.path.join(currentPath,args["output"]), args["objectiveComponent"], globalSamp, args["timeLimit"], args["gap"], loci, args["ptdm"])
#if args["globalOption"] == "mixed":
#    if args["sample"] != "all":
#        pf.strainSolver(currentPath+"/variantsAndProp/{}".format(args["sample"]), ref_strains, currentPath+"/"+args["output"], loci, args["objectiveComponent"], args["sample"], args["timeLimit"], args["gap"])
#    else:
#        pf.strainSolver(currentPath+"/variantsAndProp", ref_strains, currentPath+"/"+args["output"], loci, args["objectiveComponent"], args["sample"], args["timeLimit"], args["gap"])
#else:
#    if args["sample"] != "all":
#        dataPath = currentPath+"/variantsAndProp/{}".format(args["sample"])
#    else:
#        dataPath = currentPath+"/variantsAndProp"

    
#    print("\n~~~~~~~~~~~~~~~~~~~~~~ Solving ILP for strain prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
#    solution_dict, ilp_objective_dict, data, strains, newNameToOriName = pf.minNewStrain(dataPath, ref_strains, loci, args["sample"])
#    #    print(solution_dict)
#    print("Number of solutions from ILP: {}".format(len(solution_dict)))
#    print("\n~~~~~~~~~~~~~~~~~~~~~~ Solving LP for proportion prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
#    feasible_sol = list()
#    lp_objective_dict = dict()
#    sample_strainProp_dict = dict()
#    infeasibility = 0
#    for i in solution_dict.keys():
#        feasible = False
#        print("============= ILP Solution {0} ===============".format(i))
#        try:
#            objvalue,errObj, propObj, sampleAndStrainProp, feasible = pf.minNewStrainProp(solution_dict[i], data, strains, ref_strains, loci, newNameToOriName)
#            if feasible == False:
#                infeasibility += 1
#        except cplex.exceptions.errors.CplexSolverError as e:
#            infeasibility +=1
#            print("Infeasible!")
#            
#        if feasible == True:
#            print("Objective Value (Strain, Prop, Error): ({0}, {1}, {2})".format(ilp_objective_dict[i],propObj , errObj))
#            print("Total objective: {}".format(objvalue+ilp_objective_dict[i]))
#            feasible_sol.append(i)
#            lp_objective_dict[i] = objvalue
#            #sampleAndStrainProp is a dict with samples as keys and dataframes as values
#            sample_strainProp_dict[i] = sampleAndStrainProp
#        
#    if infeasibility == len(solution_dict):
#        print("No feasible solutions")
#    else:
#        total_obj_dict = dict()
#        min_obj = np.inf
#        for j in feasible_sol:
#            total_obj_dict[j] = ilp_objective_dict[j] + lp_objective_dict[j]
#            if total_obj_dict[j] < min_obj:
#                min_obj = total_obj_dict[j]
#        
#        min_obj_list = [i for i in total_obj_dict.keys() if total_obj_dict[i] == min_obj]
#        print("\nMinimum objective value: {}".format(min_obj))
#        print("\nSolution which minimizes both pure ILP and LP: {}\n".format(min_obj_list))
#        if len(min_obj_list) > 1:
#            print("@@@@@@@@@@@@@@@ More than 1 optimal solution for the problem @@@@@@@@@@@@@@@@@@")
#            
#        for samp in sample_strainProp_dict[min_obj_list[0]].keys():
#            output = sample_strainProp_dict[min_obj_list[0]][samp]
#            print output
#            output.to_csv("{0}/{1}_strainsAndProportions.csv".format(currentPath+"/"+args["output"], samp))
    
print("")
print("Script done.")
print("Created folder strainsAndProp which contains strains and their proportions for each sample")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
