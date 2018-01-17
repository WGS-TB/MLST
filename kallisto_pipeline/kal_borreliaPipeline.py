import os
import kal_functions as kf
import numpy as np
import itertools
import pandas as pd
import csv
import time
import argparse

start_time = time.time()

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

#variant prediction    
ap = argparse.ArgumentParser()
ap.add_argument("-d", "--data", required=True)
ap.add_argument("-c", "--numOfCores", default=4, required=False)
ap.add_argument("-f", "--objectiveComponent", default="all", required=False)
args = vars(ap.parse_args())
#variantPred_cmd = "python kal_variantAndPropPrediction.py -d {0} -c {1}".format(args["data"], args["numOfCores"])

print("")
print("******************** Solving ILP to find strains and their proportions in each sample **************************")
print("")
kf.strainSolver(currentPath+"/variantsAndProp", currentPath+"/strain_ref.txt", currentPath+"/strainsAndProp", loci, args["objectiveComponent"])

print("")
print("Script done.")
print("Created folder strainsAndProp which contains strains and their proportions for each sample")
print("Time taken : {} hr(s)".format((time.time() - start_time)/3600))
