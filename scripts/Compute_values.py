#!/usr/bin/python
from __future__ import division
from scipy.stats import binom
from collections import defaultdict
from scipy.special import comb
from pandas import *
import operator
import subprocess
import sh
import math
import random
import os
import argparse
import matplotlib.pyplot as plt
import csv
from itertools import groupby
import matplotlib.pyplot as plt
import variantILP as ilp


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

def upperfirst(x):
    return x[0].upper() + x[1:]

def precision(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(predicted))
        
def recall(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(true))

def tree():
    return defaultdict(tree)

def predictedCorrectly(predicted, true):
        return set(predicted) == set(true)

#function to generate matrix to be passed to 
def Generate_Matrix(path):
        var_list = [] #holds the variants
        read_list = [] #holds the reads
        mismatch_list = [] #holds the mismatches
        with open(path) as inf:
                for line in inf:
                        parts = line.split('\t') # split line into parts
                        if len(parts) > 1:   # if at least 2 parts/columns
                                read_list.append(parts[0]) #append reads to a list
                                var_list.append(parts[1])   #append vars to a list
                                mismatch_list.append(parts[2]) #append mismatches to list
                flag = True
        if flag is True:
                d = defaultdict(list) #dictionary holding all the variants that a read maps to
                d_1 = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
                d_2 = defaultdict(list) #dictionary holding indices for later use


                for i in range(len(read_list)):
                        num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
                        d[read_list[i]].append(var_list[i]) #append all the variants that read read_i maps to
                        d_1[read_list[i]].append(num_mismatch) #append all the mismatches that each read_i has when it maps to a variants
                        d_2[read_list[i]].append(i) #for testing purposes

                var_list = set(var_list)
                x = tree()
 
#               create a 2D dictionary that contains all the possible combinations of a read with a variant and the number of mismatches.
                for var in var_list:
                        for key in d:
                                temp = d[key]
                                if var in temp:
                                        index = temp.index(var)
                                        val = d_1[key][index]
                                        #only for perfect matches.
                                        #if val == 1:
                                        #print val
                                        x[key][var] = int(val)

                df = DataFrame(x).T.fillna(-1)
        return df

def compute_probability(n, k):
        if k > 3:
                return 0
        else:
                #b = comb(n, k, exact=False)
                x = math.pow(0.99,(n-k))
                y = math.pow(0.01,k)
                prob = x*y
        return prob

def compute_proportions(dataframe):
        prob_list = [] #a list to hold the probabilities
        for row in dataframe.itertuples(index=False):
                temp_list = list(row)
                run_sum = 0
                #compute the probability for each row in the matrix
                for i in range(len(temp_list)):
                        if temp_list[i] >= 0:
                                temp_list[i] = compute_probability(76,int(temp_list[i]))
                        else:
                                temp_list[i] = 0
                total = sum(temp_list)
                #solve for k
                temp_list = [j*(1.0/total) for j in temp_list]
                prob_list.append(temp_list)
        col_sums = [sum(k) for k in zip(*prob_list)]
        total_sum = sum(col_sums)
        prop_list = [l*(1/total_sum) for l in col_sums]
        return prop_list     

def compute_True_objective_val(dataframe):
        run_sum = 0
        for row in dataframe.itertuples(index=False):
                my_list = list(row)
                run_sum += min(i for i in my_list if i >= 0)
        run_sum += len(variants_current)
        return run_sum

def create_dictionary(keys, vals):
        my_dict = dict()
        if len(keys) == len(vals):
                for i in range(len(keys)):
                        my_dict[keys[i]] = vals[i]*100
        return my_dict 

ap = argparse.ArgumentParser()
ap.add_argument("-g", "--gene", required = True,  help="name of gene")
args = vars(ap.parse_args())
random.seed(a=1994) #set random seed
true_ratios = []
true_variants = []
total_var = []
Precision = []
Recall = []
pred_object_vals = []
true_Objective_vals = []
count = 0
for x in range(1,21): 
        k =random.randint(2,7) #generate a random integer k between 2 and 7
	#generate k random fractions that sum up to 1
        r = [random.random() for j in range(k)] 
        s = sum(r)
        r = [ i/s for i in r ]
	#run a bash sub command to sort the variants text file and return the firt k variants
        variants = (sh.head(sh.sort("variants.txt", "-R"), "-n", k))
        variants = list(variants) #convert output from runningCommand to a list
        #print variants #for testing purposes
        ratios = [] #list to store the proprotions of reads to generate
        true_prop = dict()
	gene = args["gene"]
	w = gene[0].capitalize()
	gene = gene.replace(gene[0],w)
	temp2 = ""
	file_name2 = gene+"_"+str(x)+".fas"
        variants_current = [] #list to hold the variants randomly chosen
        for i in range(len(r)):
                
                string = str(variants[i]) 
                string = string.rstrip() #remove the "\n" character that is returned by bash
                string1= string.split(">")
		#print string1
                variants_current.append(string1[1]) #append the variant to the list
                sim_name = string1[1]+"_"+str(x)+"_"
                file_name = sim_name+"reference.fa"
                val = math.ceil(r[i]*197) #compute the proportions 
                ratios.append(int(val)) 
                temp=sh.grep(string,"linear.txt","-w","-A1") #use bash to extract the variant sequence
		temp = temp.rstrip()
                temp = str(temp)
		temp2 = temp2 + temp
		#print temp2 
                cmd = "art_illumina -ss HS25 -sam -i "+file_name+" -p -l 76 -c "+str(ratios[i])+" -m 200 -s 30 -o "+sim_name #the ART command to generate the simulated data for a variant
                #write the variant sequence to a text file
		with open(file_name, "w") as text_file:
                        text_file.write(temp)
                #print sim_name #testing purposes
                #print file_name #testing purposes
                #print cmd #testing purposes
                os.system(cmd) #run the ART command
	with open(file_name2, "w") as text_file2:
		text_file2.write(temp2)
        new_cmd = "cat *_1.fq > "+gene+"_"+str(x)+"_1.fa" #append all the first of the pairs together
        new_cmd2 ="cat *_2.fq > "+gene+"_"+str(x)+"_2.fa" #append all the second of the pairs together
	new_cmd3 = "bash /home/elijah/Desktop/SRA_bowtie/Alignments/scripts/temp.sh "+ gene+"_"+str(x) + " " + gene + "_" + str(x)
	#print new_cmd3 
        os.system(new_cmd) #run the command
        os.system(new_cmd2) #run the command
	os.system(new_cmd3)
        os.system("rm "+args["gene"]+"*") #remove unneccessary files for the next iteration.
        #print ratios #for testing purposes
        X =  gene+"_"+str(x)
        X = upperfirst(X)
        #print X
        
        true_ratios.append(r)
        true_variants.append(variants_current)
        for j in range(0,len(variants_current)):
                key = variants_current[j]
                true_prop[key] = float(r[j])*100
        
        path = X+'_reads.txt'
        df = Generate_Matrix(path)
        df.rename(columns={'Unnamed: 0': 'Read'}, inplace=True)
        pred_object_val,var_predicted,reads_cov = ilp.solver(df)
        
        Precision.append(precision(var_predicted, variants_current))
        Recall.append(recall(var_predicted, variants_current))
        pred_object_vals.append(pred_object_val)
        df1 = df[var_predicted]
        df2 = df.loc[reads_cov,variants_current]
        df3 = df.loc[reads_cov,var_predicted]
        prop = compute_proportions(df3)
        pred_prop = create_dictionary(var_predicted, prop)
        val = totalVariationDist(pred_prop, true_prop)
        total_var.append(val)
        print 'True variants are:', variants_current
        print 'Predicted variants are:', var_predicted
        print 'True proportions are:', true_prop
        print 'Predicted proportions are:', pred_prop
        true_Objective_vals.append(compute_True_objective_val(df2))

        if predictedCorrectly(var_predicted, variants_current):
                count += 1

print "total variation distances are:",total_var
avg = sum(total_var)/len(total_var)   
variance = map(lambda x: (x - avg)**2, total_var)
variance = sum(variance)/len(variance)
stdev = math.sqrt(variance)
print "The average variation distance is:", avg
print "The variance is:", variance
print "standard deviation is:",stdev
print 'Precision is:', Precision
print 'Recall is:', Recall
print 'Predicted objective values are:', pred_object_vals
print 'True objective values are:', true_Objective_vals
print x
print count/x
X = []
for i in range(1,21):
        X.append(i)
plt.hist(total_var, bins='auto')
plt.xlabel('Total Variation Distance')
plt.ylabel('Frequency')
plt.show()
