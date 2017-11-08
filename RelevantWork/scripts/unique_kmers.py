from __future__ import division
from collections import defaultdict
import argparse
import subprocess
import sh
import re
#from collections import Counter
import itertools

#function to compute the kmers of a given string
def compute_kmers(string, k):
  kmers = []
  n = len(string)
  for i in range(0, n-k+1):
    kmers.append(string[i:i+k])
  kmers = set(kmers)
  return kmers

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--path", required = True, help = "Path to data table containing the list of variants")
ap.add_argument("-r", "--reads", required = True,  help="Path to file containing the reads")
ap.add_argument("-k", "--kmer_length", required = True,  help="The length of desired kmers")
args = vars(ap.parse_args())

var_list = [] #list to hold variants
read_list = [] #list to hold reads
kmers_dict = defaultdict(set)
my_dict = {}
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
  my_dict[var] = compute_kmers(sequence, k)
  all_kmers.append(compute_kmers(sequence, k))

  #print (kmers_list)
  #add the kmers to a dictionary
print "Total length of variant dictionary: ", len(all_kmers)
#print my_dict.keys()

flattened = [item for sublist in all_kmers for item in sublist]
print "Total number of kmers:", len(flattened)
unique = set(flattened)
print "Total number of unique kmers:", len(unique)

for kmer in unique:
  for key in my_dict:
    if kmer in my_dict[key]:
      kmers_dict[kmer].update([key])
      print key

unique_dict = {}
count = 0
print "kmer,", "unique variant"
for kmer in kmers_dict:
  if len(kmers_dict[kmer]) == 1:
    count += 1
    print kmer, kmers_dict[kmer]
    unique_dict[kmer] = kmers_dict[kmer]
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
	print kmer, prop_dict[kmer]
