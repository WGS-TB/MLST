#!/usr/bin/python
from __future__ import division
import argparse
import subprocess
import sh
import os
import re


#function to compute the kmers of a given string
def compute_kmers(string, k):
  kmers = []
  n = len(string)
  for i in range(0, n-k+1):
    kmers.append(string[i:i+k])
  kmers = set(kmers)
  return kmers

# construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--path", required = True, help = "Path to data table containing the list of variants")
ap.add_argument("-r", "--reads", required = True,  help="Path to file containing the reads")
ap.add_argument("-k", "--kmer_length", required = True,  help="The length of desired kmers")
ap.add_argument("-s", "--sim_var", required = True, help = "Path to data table containing the list of simulated variants")
args = vars(ap.parse_args())
  
#my_list = compute_kmers(string, 15)
#print my_list

var_dict = {} #dict to hold the variants and their kmers
var_list = [] #list of variants
sim_var_list = [] #list to hold the simulated variants
var_nums = {} #number of kmers in variants
var_prop = {} #kmers proportions
var_not_present = {} #dictionary of variants not present

#append the variants to a list
with open(args["path"]) as inf:
    for line in inf:
        var_list.append(line)


#append the simulated variants to a list
with open(args["sim_var"]) as inf3:
    for line_3 in inf3:
        sim_var_list.append(line_3.rstrip())
print sim_var_list

reads = []

with open(args["reads"]) as inf2:
  for line_2 in inf2:
    parts = line_2.split("\t")
    reads.append(parts[3])

k = int(args["kmer_length"])
read_kmers = []
#compute the set of kmers for all reads
for read in reads:
    read = read.rstrip()
    #compute the kmers
    read_kmers.append(compute_kmers(read, k))
print type(read_kmers)
flattened = [val for sublist in read_kmers for val in sublist]
print len(flattened)
#remove duplicates from the list of kmers
unique_kmers = set(flattened)
print len(unique_kmers)

#generate kmers for the variants 
for var in var_list:
  #print var
  var = var.rstrip()
  #read_list = []
  read_list = (sh.grep(var,"linear.txt","-w","-A1")) #get the read
  read_list = str(read_list)
  #process the read
  read_list = re.findall('\d*\D+',read_list)
  temp = read_list[1] 
  temp = temp.rstrip()
  #remove new lines
  temp = temp.split("\n")
  #extract the sequence of interest
  sequence = temp[1]
  #compute the kmers for a given read
  var_kmers = compute_kmers(sequence, k)
  #print (kmers_list)
  #add the kmers to a dictionary
  if var in var_dict:
    var_dict[var].append(var_kmers)
  else:
    var_dict[var] = var_kmers


for kmer in var_dict:
  kmers = var_dict[kmer]
  #print len(kmers)
  common = set.intersection(kmers, unique_kmers)
  prop = len(common)/len(kmers)
  smallest = 100
  #
  smallest = min(smallest,prop)
  var_nums[kmer] = len(common)
  var_prop[kmer] = round(prop, 4)
#print var_nums
#print var_prop
print smallest
with open("output_k_35.txt", 'w') as f:
  for key, value in var_prop.items():
    f.write('%s:%s\n' % (key, value))

#get the variants not used for the simulated reads
for kmer in var_prop:
  if var_prop[kmer] >= smallest:
    var_not_present[kmer] = var_prop[kmer]
print len(var_not_present)
print len(var_prop)

with open("output_k_35_S.txt", 'w') as g:
  for key, value in var_not_present.items():
    g.write('%s\n' % (key))