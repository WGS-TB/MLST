#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 14:54:38 2018

@author: Elijah Willie

The below code generates a fasta file containing the sequences of each strain within the MLST schema
"""

#import the required libraries
from __future__ import division
from collections import defaultdict
import pandas as pd
import os
import random
import sh
import csv
import sys
import numpy as np

#get the currenct working directory 
curr_dir = os.getcwd()

#get the folder containing the variant sequences
variant_files = os.path.join(curr_dir, 'Variant_files')

#read in the MLST schema
loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
schema = pd.read_csv(os.path.join(curr_dir, 'strain_ref.txt'), sep = '\t', usecols=range(1,len(loci)+1))

for name in loci:
        schema["%s" %name] = name + "_" + schema["%s" %name].astype(str)

#now for each strain, extract all the allele sequences and append it to a fasta file
#change into the appropriate directory
os.chdir(variant_files)
my_str = open('Names.csv','w')
count = 1
for index in range(schema.shape[0]):
    count = count + 1
    #open a file to append the sequences onto
    sequence_file = open('Strain_{}_Genome.fas'.format(index), 'w')
    my_str.write('Strain_{}_Genome.fas '.format(index))
    #extract the current strain
    row = schema.iloc[index,:].values.tolist()
    #iterate over all the alleles in the current strain
    for allele in row:
        #extract the allele sequence
        allele_sequence = sh.grep(allele,"Genome_linear.fasta","-w","-A1")
        #remove the newline (\n) character
        allele_sequence = allele_sequence.rstrip()
        #convert the sequence back to a string
        allele_sequence = str(allele_sequence)
        #write this sequence to the current genome file
        sequence_file.write('{0}{1}'.format(allele_sequence, '\n'))
    
    #close the sequence file
    sequence_file.close()
    if count == 100:
        break
my_str.close()

