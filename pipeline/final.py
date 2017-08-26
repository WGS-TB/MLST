import os
import getVariantAndProp as gvp
import getStrAndProp as gsp
import numpy as np
import itertools
import pandas as pd
import csv
import time
import argparse

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

if not os.path.exists("strainsAndProp"):
    os.mkdir("strainsAndProp")

print("")
print("******************** Solving ILP to find strains and their proportions in each sample **************************")
print("")
#if args["sample"] != "all":
#    gsp.strainSolver(currentPath+"/variantsAndProp/{}".format(args["sample"]), currentPath+"/strain_ref.txt", currentPath+"/strainsAndProp", loci, "all", args["sample"])
#else:
#    gsp.strainSolver(currentPath+"/variantsAndProp", currentPath+"/strain_ref.txt", currentPath+"/strainsAndProp", loci, "all", args["sample"])
if args["sample"] != "all":
    print("~~~~~~~~~~~~~~~~~~~~~~ Solving ILP for strain prediction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    solution_dict, data, strains = gsp.minNewStrain(currentPath+"/variantsAndProp/{}".format(args["sample"]), currentPath+"/strain_ref.txt", loci, args["sample"])
    print("Number of solutions from ILP: {}".format(len(solution_dict)))
    print("~~~~~~~~~~~~~~~~~~~~~~~~ Solving LP for proportion prediction ~~~~~~~~~~~~~~~~~~~~~~~~~")
    for i in solution_dict.keys():
        gsp.minNewStrainProp(solution_dict[i], data, strains, ref_strains, loci)
else:
    solution_dict, data, strains = gsp.minNewStrain(currentPath+"/variantsAndProp".format(args["sample"]), currentPath+"/strain_ref.txt", loci, args["sample"])
    
    for i in solution_dict.keys():
        gsp.minNewStrainProp(solution_dict[i], data, strains, ref_strains, loci)

print("")
print("Script done.")
print("Created folder variantsAndProp which contains variants identified and their proportions for each sample")
print("Created folder strainsAndProp which contains strains and their proportions for each sample")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
