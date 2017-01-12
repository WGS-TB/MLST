#!/usr/bin/python
from scipy.stats import binom
from scipy.special import comb
import math
import argparse
import csv

#function to compute the binomial probability of error in a read
#where n is the read length and k is the number of mismatches
def compute_probability(n, k):
	b = comb(n, k, exact=False)
	x = math.pow(0.99,(n-k))
	y = math.pow(0.01,k)
	prob = b*x*y
	return prob

# construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--path", required = True, help = "Path to data table")
ap.add_argument("-g", "--gene", required = True,  help="name of gene")
args = vars(ap.parse_args())

var_list = [] #holds the variants
read_list = [] #holds the reads
mismatch_list = [] #holds the mismatches
read_dict = dict() #dictionary to hold the reads and their total probabilities
prop_dict = dict() #dictionary to hold the variants and their proportions
var_dict = dict() #dictionary to hold just the variants
k_dict = dict() #dictionary to hold the k-values for the reads

#open table for reading
with open(args["path"]) as inf:
    for line in inf:
        parts = line.split('\t') # split line into parts
        if len(parts) > 1:   # if at least 2 parts/columns
	    read_list.append(parts[0]) #append reads to a list
            var_list.append(parts[1])   #append vars to a list
            mismatch_list.append(parts[2]) #append mismatches to list

#compute variant dictionary
count = 0 #don't worry about the count. It's just a place holder for the dictionary
for i in xrange(0, len(var_list)):
	key = var_list[i]
	if not (var_dict.has_key(key)):
		var_dict[key] = count
		count += 1
print 'length variant dictionary is:', len(var_dict)

#compute P[read_j|var_i] for all i
for i in xrange(0, len(read_list)): #iterate over all the reads
	key = read_list[i]
	num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
	prob = compute_probability(76, num_mismatch) #compute the binomial probability
	#check if this variant probability has been added for a read
	if not (read_dict.has_key(key)): #if read is not in the dictionary, add it, and also add its probability
		read_dict[key] = prob
	else: #if it's already in the dictionary, just update its probability
		read_dict[key] += prob
print 'length read dictionary is:', len(read_dict)

#compute the k_values by dividing 1 by the sum of probabilities of the the variants that each read maps to
for key in read_dict:
	val = read_dict.get(key)
	k = 1/(val) #compute the k_j value
	k_dict[key] = k #store in in the dictionary
print 'length k_j dictionary is:', len(k_dict)

#compute the sum of probabilities of all reads that map to a variant
#this can be computed for var_i as follows:
#prop = k_1*P[read_1|var_i] + k_2*P[read_2|var_i] +....+ k_j*P[read_j|var_i]
for key in var_dict:
	for j in xrange(0, len(var_list)):
		temp = read_list[j]
		k_j = k_dict.get(temp)
		if key == var_list[j]:
			num_mismatch = mismatch_list[j].count('>') #count the number of mismatches for each read
        		prob = compute_probability(76, num_mismatch) #compute the binomial probability
			if not (prop_dict.has_key(key)): #if key not in dictionary, insert it and then add its first k_j*P[read_j|var_i]
				prop_dict[key] = k_j*prob 
			else: #if variant is already in the dictionary, just add k_j*P[read_j|var_i] to the growing sum
				prop_dict[key] += (k_j*prob)
#we should now have a dictionary (prop_dict) where each (key, value) pair is a variant and the sum of all the probabilites of reads
#that maps to it. Where each term in the sum is of the form k_j*P[read_j|var_i]
print 'length proportion dictionary is:', len(prop_dict)

#now sum over all (key, value) pair in prop_dict to compute the total sum
total_sum = sum(prop_dict.values()) #sum all probability values in the list
print 'total sum is:', total_sum

#now to calculate the proportions of each variant
for key in prop_dict:
	var_prop = (prop_dict.get(key)/total_sum)*100 #the proportion of a variant
	prop_dict[key] = var_prop #update the key value to have the proportions
#write thes proportions to a table
w = csv.writer(open(args["gene"]+'_proportions.csv', "w"))
for key, val in prop_dict.items():
    w.writerow([key, val])	
