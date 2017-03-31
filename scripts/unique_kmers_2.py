from __future__ import division
from collections import defaultdict
import argparse
import subprocess
import sh
import re
from collections import Counter
import itertools

#function to compute the kmers of a given string
def compute_kmers(string, k):
  kmers = []
  n = len(string)
  for i in range(0, n-k+1):
    kmers.append(string[i:i+k])
  #kmers = set(kmers)
  return kmers

def generate_pairs(k_list):
  return [(k_list[x1], k_list[x2]) for x1 in range(len(k_list)) for x2 in range(x1+1, len(k_list))]

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--path", required = True, help = "Path to data table containing the list of variants")
ap.add_argument("-r", "--reads", required = True,  help="Path to file containing the reads")
ap.add_argument("-k", "--kmer_length", required = True,  help="The length of desired kmers")
args = vars(ap.parse_args())

var_list = [] #list to hold variants
read_list = [] #list to hold reads
kmers_dict = defaultdict(set)
my_dict_pairs = defaultdict(set)
my_dict = {}
all_pairs = []
all_kmers = []


#append the variants to a list
with open(args["path"]) as inf:
    for line in inf:
        var_list.append(line)

with open(args["reads"]) as inf2:
  for line_2 in inf2:
    parts = line_2.split("\t")
    read_list.append(parts[3])

k = int(args["kmer_length"])

#generate kmers for the variants 
for var in var_list:
  #print var
  var = var.rstrip()
  #read_list = []
  my_seq = (sh.grep(var,"linear.txt","-w","-A1")) #get the read
  my_seq = str(my_seq)
  #process the read
  my_seq = re.findall('\d*\D+',my_seq)
  temp = my_seq[1] 
  temp = temp.rstrip()
  #remove new lines
  temp = temp.split("\n")
  #extract the sequence of interest
  sequence = temp[1]
  #compute the kmers for a given read
  #my_dict[var] = compute_kmers(sequence, k)
  #all_kmers.append(compute_kmers(sequence, k))
  #compute the unique pairs
  current = compute_kmers(sequence, k)
  y = set(generate_pairs(current))
  my_dict_pairs[var] = y
  all_pairs.append(generate_pairs(current))


flattened = [item for sublist in all_pairs for item in sublist]
unique = set(flattened)
  #print (kmers_list)
  #add the kmers to a dictionary
#print "Total length of variant dictionary: ", len(var_list)
#print my_dict.keys()
#print "Total number of pairs of kmers:", len(flattened)

#print "Total number of unique pairs of kmers:", len(unique)

unique_2 = {}
currentSet = set()
total = 0
for var_1 in var_list:
  var_1= var_1.rstrip()
  currentSet = my_dict_pairs[var_1]
  for var_2 in var_list:
    var_2 = var_2.rstrip()
    if var_2 != var_1:
      currentSet.difference_update(my_dict_pairs[var_2])
  total += len(currentSet)
  unique_2[var_1] = currentSet
print total
len(unique_2)

count_dict = {}

print "variant", "Number of unique pairs"
for var in unique_2:
  if len(unique_2[var]) > 0:
    print var, len(unique_2[var])
    count_dict[var] = len(unique_2[var])

prop_dict = {}
print "variant", "proportion of unique pairs"
for var in count_dict:
  prop_dict[var] = round((count_dict[var]/total)*100, 4)
  print var, round((count_dict[var]/total)*100, 4)


'''for kmer_pair in unique:
  for key in my_dict_pairs:
    if kmer_pair in my_dict_pairs[key]:
      kmers_dict[kmer_pair].update([key])

unique_dict = {}
count = 0
print "kmer_pair,", "unique variant"
for kmer_pair in kmers_dict:
  if len(kmers_dict[kmer_pair]) == 1:
    count += 1
    print kmer_pair, kmers_dict[kmer_pair]
    unique_dict[kmer_pair] = kmers_dict[kmer_pair]
print "Total number kmers that maps to one variant:", count

prop_dict = {}
count_dict = {}
#compute the counts from proportions
for kmer in unique_dict:
  x = list(unique_dict[kmer])
  if not (count_dict.has_key(x[0])):
    count_dict[x[0]] = 1
  else:
    count_dict[x[0]] += 1
print "Total number of variants with at least one unique kmer:", len(count_dict)
print "variant,", "counts of unique kmers"
for kmer in count_dict:
  print kmer, count_dict[kmer]

#compute the proportions
for kmer in count_dict:
  prop_dict[kmer] = (count_dict[kmer]/count)*100

print "variant,", "counts of unique kmers"
for kmer in prop_dict:
  print kmer, prop_dict[kmer]'''
