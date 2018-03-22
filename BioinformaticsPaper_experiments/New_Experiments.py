#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 14:20:23 2018

@author: elijah
"""

from __future__ import division
from collections import defaultdict
import pandas as pd
import csv
import numpy as np
import random
import os
import matplotlib.pyplot as plt
import pipeline_functions as pf
import argparse
import itertools
import time
import numpy as np, numpy.random
import variantILP as varSolver
import sh
plt.style.use('ggplot')

#Set seed for reproducibility

seed = 1995
random.seed(seed)

#set up genes edit distance matrices

#clpA
clpA_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_clpA.csv', sep=",")
#clpX
clpX_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_clpX.csv', sep=",")
#nifS
nifS_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_nifS.csv', sep=",")
#pepX
pepX_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_pepX.csv', sep=",")
#pyrG
pyrG_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_pyrG.csv', sep=",")
#recG
recG_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_recG.csv', sep=",")
#rplB
rplB_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_rplB.csv', sep=",")
#uvrA
uvrA_df = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/editDistanceMatrix_uvrA.csv', sep=",")

#read in the reference strain database
loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
genes_df = [clpA_df,clpX_df,nifS_df,pepX_df,pyrG_df,recG_df,rplB_df,uvrA_df]
reference = pd.read_csv('/home/elijah/Documents/Borellia/New_Experiments/New_strain_ref.txt',sep="\t",usecols=range(1,len(loci)+1))
#parse the database
for name in loci:
    reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)

#add all the strains to a list    
all_strains = []
for index in range(reference.shape[0]):
    all_strains.append(reference.iloc[index,:].values.tolist())

#check if a strain is already in the database    
def check_strain(strain,strains):
    if strain not in strains:
        return True
    else:
        return False
        
def Dict_to_csv(gene_dict, filename):
    '''
        A function that takes in a dictionary and writes it to a csv file with the given name

        gene_dict: The dictionary to be written

        filename: The name of the file
    '''
    with open(filename+'_proportions.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in gene_dict.items():
            writer.writerow([key, value])
        
#function that computes the allele that is within max_editDistance of selected allele
def Compute_Allele(A, max_editDistance):
    allele_new = A.split('_')[0]
    idx = loci.index(allele_new)
    allele_df = genes_df[idx]
    allele_list = allele_df.loc[allele_df['level_0'] == A]
    allele_list = allele_list
    allele_dict = allele_list.to_dict()
    allele_dict.pop(A)
    max_dist_alleles = []
    for allele, dist in allele_dict.iteritems():    # for name, age in list.items():  (for Python 3.x)
       if list(dist.items()[0])[1] <= max_editDistance and allele != A:
           max_dist_alleles.append(allele)
    if len(max_dist_alleles) > 1:
        return random.choice(max_dist_alleles)
    if len(max_dist_alleles) == 0:
        return min(allele_dict, key=allele_dict.get)
    return max_dist_alleles[0]


def Mutate_strain(S1, editDist, num_mut):
    S2 = [None]*8
    for j in range(len(S1)):
        S2[j] = S1[j]
    #radomly select an allele to mutate
    for i in range(num_mut):    
        A = random.choice(S1)
        #find closest allele with max edit distance
        L = Compute_Allele(A, editDist)
        idx = S1.index(A)
        #replace this allele to create a new strain S2
        S2[idx] = L
    #print 'strain one is:', S1
    #print 'strain two is:', S2
    return S2

def Recombine_strains(strain1, strain2):
    result = []
    for i in range(1,9):
        pos = np.random.choice(np.arange(1, 3), p=[0.5,0.5])
        if pos == 1:
            result.append(strain1[i-1])
        else:
            result.append(strain2[i-1])
    return result
    


def Dict_to_csv(gene_dict, filename):
    '''
        A function that takes in a dictionary and writes it to a csv file with the given name

        gene_dict: The dictionary to be written

        filename: The name of the file
    '''
    with open(filename+'_proportions.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in gene_dict.items():
            writer.writerow([key, value])

def writeReadTable():
    readOutFile = open("all_Strains_pairedNoHeader.sam")
    writefile = open("all_Strains_reads.txt", "w")
    for line in readOutFile:
        fields = line.strip("\t").split()
        read = fields[0]
        allele = fields[2]
        quality = fields[10]
#        mm = [i for i in fields if i.startswith("XM:i:")][0]  #bowtie2
        mm = [i for i in fields if i.startswith("NM:i:")][0]   #bowtie
        mm_pos = [j for j in fields if j.startswith("MD:Z:")][0]
        
        writefile.write(read + "\t" + allele + "\t" + quality + "\t" + mm + "\t" + mm_pos + '\n')
        
    readOutFile.close()
    writefile.close()
            
def Compute_Variant_proportions(strain,strain_prop):
    '''
        A function that takes in a strain, its proportion and computes the proportions
        of its alleles

        strain: The current strain under scrutiny

        strain_prop: The proportions to be assigned to its alleles
    '''
    variant_proportions = defaultdict(list)
    for gene in strain:
        variant_proportions[gene] = strain_prop#*100
    #print sum(variant_proportions.values())
    return variant_proportions

def Dict_to_csv(gene_dict, filename):
    '''
        A function that takes in a dictionary and writes it to a csv file with the given name

        gene_dict: The dictionary to be written

        filename: The name of the file
    '''
    with open(filename+'_proportions.csv', 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in gene_dict.items():
            writer.writerow([key, value])


            
    
def Generate_reads(strains, seed, iteration, editDist):
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    nums = [random.uniform(0,1) for x in range(0,len(strains))]
    Sum = reduce(lambda x,y: x+y, nums)
    proportions = [x/Sum for x in nums]
    for i in range(len(strains)):
        sequence_file = open('editDist_{}_iteration_{}_Strain_{}_sequences.fas'.format(editDist,iteration,i), 'w')
        strain = strains[i]
        output_file_name = 'editDist_{}_iteration_{}_strain_{}_reads_'.format(editDist,iteration,i)
        for j in range(len(strain)):
            variant_sequence = sh.grep(str(strain[j]),"Variant_files/{}_linear.txt".format(genes[j]),"-w","-A1") #use bash to extract the variant sequence
            variant_sequence = variant_sequence.rstrip() #remove the "\n" character that is returned by bash
            variant_sequence = str(variant_sequence)
            #total_variants_sequence += variant_sequence
            sequence_file.write('{0}{1}'.format(variant_sequence, '\n'))
        #run art to generate reads with 100X coverage
        #store
        sequence_file.close()
        #Set the ART command, I have included a random seed for reproducibility, and a coverage parameter
        ART_command = "art_illumina -qL 33 -qs 10 -qs2 15 -k 3 -rs {0} -q -ss HS25 -sam -i {1} -p -l 76 -f {2} -m 200 -s 10 -o {3}".format(seed, 'editDist_{}_iteration_{}_Strain_{}_sequences.fas >/dev/null 2>&1'.format(editDist,iteration,i),proportions[i]*100, output_file_name)
        os.system(ART_command)
    appendFirst_cmd = "cat *_1.fq > editDist_{}_iteration_{}_all_Strains_1.fa".format(editDist,iteration) #append all the first of the pairs together
    appendSecond_cmd ="cat *_2.fq > editDist_{}_iteration_{}_all_Strains_2.fa".format(editDist,iteration) #append all the second of the pairs together
    os.system(appendFirst_cmd)
    os.system(appendSecond_cmd)
    os.system('rm *aln*')
    os.system('rm *seq*')
    os.system('rm *sam*')
    #os.system('rm *fa*')
    os.system('rm *fq*')
    return proportions    

    
def Compute_ADP_Prec_and_rec(true, predicted):
    '''
    A function that computes the precision and recall afer a simulation
    
    true: A list that contains the true alleles
    
    Predicted: A list that contians the predicted alleles from the ADP
    
    Returns a precision and recall value
    '''
    tvd = 0
    count = 0
    print ('The True Alleles  and their proportions are are {}\n'.format(true))
    print ('The Predicted Alleles and their proportions are {}\n'.format(predicted))
    
    for key in true.keys():
        if key in predicted.keys():
            count = count + 1
            tvd = tvd + (true[key] - predicted[key])
        else:
            tvd = tvd + true[key]
    for key in predicted.keys():
        if key not in true.keys():
            tvd = tvd + predicted[key]
    prec = count/len(true.keys())
    rec = count/len(predicted.keys())
    #print len(true), len(predicted)
    #print 'Count is:', count
    return prec, rec, tvd/16.0
    
    
    
def Compute_Prec_and_rec(strain_dict, strain_df):
    '''
        A function to compute the precision and recall after a simulation.

        strain_dict: A dictionary containing the true strains

        strain_df: A dataframe containing the strains predicted by the ILP

    '''
    count = 0
    total_var_dist = 0
    predicted_dict = defaultdict(list)
    #convert the dataframe into a dictionary
    for i in range(strain_df.shape[0]):
        row = strain_df.iloc[i].tolist()
        temp_key = tuple(row[0:8])
        predicted_dict[temp_key] = float(row[10])
    print ('The True strains and their proportions are {}\n'.format(strain_dict))
    print ('The Predicted strains and their proportions are\n {}'.format(predicted_dict))
    #Compute the total variation distance
    for key in predicted_dict.keys():
        if key in strain_dict.keys():
            total_var_dist += abs(round(strain_dict[key],3)-predicted_dict[key])
            count += 1
        else:
            total_var_dist += predicted_dict[key]
    for key in strain_dict.keys():
	if key not in predicted_dict.keys():
	    total_var_dist += strain_dict[key]
    #compute the precision
    precision = count/strain_df.shape[0]
    #compute the recall
    recall = count/len(strain_dict)

    return precision, recall, total_var_dist/2.0
    
    
def Write_Proportions(strains, proportions, iteration, editDist, Type):
    #create a strain and its proportions dict
    strain_prop_dict = defaultdict(list)
    for i in range(len(strains)):
        #print tuple(strains[i])
        strain_prop_dict[tuple(strains[i])] = proportions[i]
    #print len(proportions)
    clpA = defaultdict(list)
    clpX = defaultdict(list)
    nifS = defaultdict(list)
    pepX = defaultdict(list)
    pyrG = defaultdict(list)
    recG = defaultdict(list)
    rplB = defaultdict(list)
    uvrA = defaultdict(list)
    
    #Add the proportions to a dictionary
    for strain_key in strain_prop_dict.keys():
            #define a dictionary to hold the current strain and its proportions
            gene_dict = defaultdict(list)
            gene_dict[strain_key] = strain_prop_dict[strain_key]
            for gene_key in gene_dict:
                #for clpA
                if 'clpA' in gene_key[0]:
                    if gene_key[0] not in clpA:
                        clpA[gene_key[0]] = gene_dict[gene_key]
                    else:
                        clpA[gene_key[0]] += gene_dict[gene_key]
                #print clpA
                #for clpX
                if 'clpX' in gene_key[1]:
                    if gene_key[1] not in clpX:
                        clpX[gene_key[1]] = gene_dict[gene_key]
                    else:
                        clpX[gene_key[1]] += gene_dict[gene_key]
                #for nifS
                if 'nifS' in gene_key[2]:
                    if gene_key[2] not in nifS:
                        nifS[gene_key[2]] = gene_dict[gene_key]
                    else:
                        nifS[gene_key[2]] += gene_dict[gene_key]
                #for pepX
                if 'pepX' in gene_key[3]:
                    if gene_key[3] not in pepX:
                        pepX[gene_key[3]] = gene_dict[gene_key]
                    else:
                        pepX[gene_key[3]] += gene_dict[gene_key]
                #for pyrG
                if 'pyrG' in gene_key[4]:
                    if gene_key[4] not in pyrG:
                        pyrG[gene_key[4]] = gene_dict[gene_key]
                    else:
                        pyrG[gene_key[4]] += gene_dict[gene_key]
                #for recG
                if 'recG' in gene_key[5]:
                    if gene_key[5] not in recG:
                        recG[gene_key[5]] = gene_dict[gene_key]
                    else:
                        recG[gene_key[5]] += gene_dict[gene_key]
                #for rplB
                if 'rplB' in gene_key[6]:
                    if gene_key[6] not in rplB:
                        rplB[gene_key[6]] = gene_dict[gene_key]
                    else:
                        rplB[gene_key[6]] += gene_dict[gene_key]
                #for uvrA
                if 'uvrA' in gene_key[7]:
                    if gene_key[7] not in uvrA:
                        uvrA[gene_key[7]] = gene_dict[gene_key]
                    else:
                        uvrA[gene_key[7]] += gene_dict[gene_key]
    #check to see if the sums are equal to 1
    '''
    print sum(clpA.values())
    print sum(clpX.values())
    print sum(nifS.values())
    print sum(pepX.values())
    print sum(pyrG.values())
    print sum(recG.values())
    print sum(rplB.values())
    print sum(uvrA.values())
    '''
    #write proportions to csv files
    dict_list = [clpA,clpX,nifS,pepX,pyrG,recG,rplB,uvrA]
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    #write the proportions to a csv file
    if not os.path.exists('{}_editDist_{}_iteration_{}_SDP_Samples'.format(Type,editDist,iteration)):
        os.mkdir('{}_editDist_{}_iteration_{}_SDP_Samples'.format(Type,editDist,iteration))
    if not os.path.exists('{}_editDist_{}_iteration_{}_SDP_output'.format(Type,editDist,iteration)):
        os.mkdir('{}_editDist_{}_iteration_{}_SDP_output'.format(Type,editDist,iteration))
    os.chdir('{}_editDist_{}_iteration_{}_SDP_Samples'.format(Type,editDist,iteration))
    for i in range(len(dict_list)):
        locus = dict_list[i]
        Dict_to_csv(locus,genes[i])
    
    os.chdir('../')
    return strain_prop_dict

    
    
    
    
'''
The first evolutionary model
EvoMod1 (Result in 2 strains)
1) Pick a strain at random from database, say S1. Let S2=S1, where S2 is the second strain.
2) For i in range(num_mut):
	- Pick a random gene g to mutate. Denote allele at gene g for S2 as A. 
	- Find an allele which is within max_editDist from A, say A'
	- Replace A with A' in S2
3) Now we have two strains (S1 and S2)
4) Assign random proportions to them
5) Input into SDP, compute Precision, TVD, Recall and number of optimal solutions from SDP

'''
def EvoMod1(reference, editDist, num_mut, iteration, Type, Etype):
    strainRef = '/home/elijah/Documents/Borellia/New_Experiments/New_strain_ref.txt'
    samplesDir = '/home/elijah/Documents/Borellia/New_Experiments/Type_1_editDist_{}_iteration_{}_SDP_Samples'.format(editDist,iteration)
    outputDir = '/home/elijah/Documents/Borellia/New_Experiments/Type_1_editDist_{}_iteration_{}_SDP_output'.format(editDist,iteration)
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    strains = []
    test = []
    #select random strain
    S1_index =  random.randint(1,696)
    S1 = reference.iloc[S1_index,:].values.tolist()
    strains.append(S1)
    #Make sure the strain is not in the database already
    flag = True
    S2 = Mutate_strain(S1, editDist, num_mut)
    while flag:
        if check_strain(S2,all_strains):
            strains.append(S2)
            flag = False
        else:
            S1 = Mutate_strain(S1, editDist, num_mut)
            flag = True
    #get the true alleles
    true_alleles = set([item for sublist in strains for item in sublist])
    true_all = defaultdict(int)
    
    #Generate the reads
    proportions = Generate_reads(strains, seed, iteration, editDist)
    for strain in strains:
        strain = list(strain)
        for gene in strain:
            if gene not in true_all:
                true_all[gene] = proportions[strains.index(strain)]
            else:
                true_all[gene] = true_all[gene] + proportions[strains.index(strain)]
    #if we are running experiment 1
    
    #if we are running the third experiment    
    if Etype == 1:
        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        strain_df = pf.strainSolver(samplesDir,strainRef,outputDir,genes,'noProp','s',10800,5)
        #compute the precision and recall
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        return prec, rec, tvd
    #if we are running the fourth experiment    
    elif Etype == 2:
        
        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~') 
        cwd = os.getcwd()
        #run the ADP algorithm
        ADP_dict, ADP_alleles, ADP_prop = run_ADP(editDist, iteration)
        strain_df = pf.strainSolver(cwd,strainRef,outputDir,genes,'noProp','s',10800,5)
        SDP_list = []
        ADP_alleles = set(item for sublist in ADP_alleles for item in sublist)
        #compute the precision and recall for the SDP output
        ADP_pred, ADP_recall, ADP_tvd = Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        return prec, rec, tvd, ADP_pred, ADP_recall, ADP_tvd
        
'''
EvoMod2 (Result in 4 strains, 2 new 2 existing)
1) Pick 2 strains at random from database
2) Repeat EvoMod1 on each of these strains, make sure that the 2 mutated strains are novel strains
3) Now we have 4 strains
4) Same stuff as above
'''

def EvoMod2(reference, editDist, num_mut, iteration, Type, Etype):
    strainRef = '/home/elijah/Documents/Borellia/New_Experiments/New_strain_ref.txt'
    samplesDir = '/home/elijah/Documents/Borellia/New_Experiments/Type_2_editDist_{}_iteration_{}_SDP_Samples'.format(editDist,iteration)
    outputDir = '/home/elijah/Documents/Borellia/New_Experiments/Type_2_editDist_{}_iteration_{}_SDP_output'.format(editDist,iteration)
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    strains = []
    #get the indices of the strains to use
    indices = random.sample(range(1, 696), 2)
    
    for index in indices:
        S1 = reference.iloc[index,:].values.tolist()
        strains.append(S1)
        flag = True
        S2 = Mutate_strain(S1, editDist, num_mut)
        while flag:
            if check_strain(S2,all_strains):
                strains.append(S2)
                flag = False
            else:
                S1 = Mutate_strain(S1, editDist, num_mut)
                flag = True
    true_all = defaultdict(int)
    
    #Generate the reads
    proportions = Generate_reads(strains, seed, iteration, editDist)
    for strain in strains:
        strain = list(strain)
        for gene in strain:
            if gene not in true_all:
                true_all[gene] = proportions[strains.index(strain)]
            else:
                true_all[gene] = true_all[gene] + proportions[strains.index(strain)]
    #print true_all
    #if we are running experiment 1
    
    #if we are running the third experiment    
    if Etype == 1:
        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        strain_df = pf.strainSolver(samplesDir,strainRef,outputDir,genes,'noProp','s',10800,5)
        #compute the precision and recall
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        return prec, rec, tvd
    #if we are running the fourth experiment    
    elif Etype == 2:
        
        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~') 
        cwd = os.getcwd()
        #run the ADP algorithm
        ADP_dict, ADP_alleles, ADP_prop = run_ADP(editDist, iteration)
        strain_df = pf.strainSolver(cwd,strainRef,outputDir,genes,'noProp','s',10800,5)
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        ADP_alleles = set(item for sublist in ADP_alleles for item in sublist)
        #compute the precision and recall for the SDP output
        ADP_pred, ADP_recall, ADP_tvd = Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        return prec, rec, tvd, ADP_pred, ADP_recall, ADP_tvd
    
    
'''
EvoMod3(Result in 5 strains, 3 new 2 existing)
1) Do EvoMod2
2) Choose any 2 strains to recombine, to produce a 3rd new strain
3) If we can afford to do, 2 recombination events. We could try that
4) ...
'''    
    
def EvoMod3(reference, editDist, num_mut, iteration, Type, Etype):
    #do EveMod2
    strainRef = '/home/elijah/Documents/Borellia/New_Experiments/New_strain_ref.txt'
    samplesDir = '/home/elijah/Documents/Borellia/New_Experiments/Type_3_editDist_{}_iteration_{}_SDP_Samples'.format(editDist,iteration)
    outputDir = '/home/elijah/Documents/Borellia/New_Experiments/Type_3_editDist_{}_iteration_{}_SDP_output'.format(editDist,iteration)
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    strains = []
    #get the indices of the strains to use
    indices = random.sample(range(1, 696), 2)
    #Mutate strains
    for index in indices:
        S1 = reference.iloc[index,:].values.tolist()
        strains.append(S1)
        flag = True
        S2 = Mutate_strain(S1, editDist, num_mut)
        while flag:
            if check_strain(S2,all_strains):
                strains.append(S2)
                flag = False
            else:
                S1 = Mutate_strain(S1, editDist, num_mut)
                flag = True
        
    #recombine strains
    S3 = Recombine_strains(reference.iloc[indices[0]].tolist(), reference.iloc[indices[1]].tolist())
    flag = True
    while flag:
        if check_strain(S3,all_strains):
            strains.append(S3)
            flag = False
        else:
            S3 = Mutate_strain(S1, editDist, num_mut)
            flag = True
    true_all = defaultdict(int)
    
    #Generate the reads
    proportions = Generate_reads(strains, seed, iteration, editDist)
    for strain in strains:
        strain = list(strain)
        for gene in strain:
            if gene not in true_all:
                true_all[gene] = proportions[strains.index(strain)]
            else:
                true_all[gene] = true_all[gene] + proportions[strains.index(strain)]
    
    #if we are running experiment 1
    
    #if we are running the third experiment    
    if Etype == 1:
        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        strain_df = pf.strainSolver(samplesDir,strainRef,outputDir,genes,'noProp','s',10800,5)
        #compute the precision and recall
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        return prec, rec, tvd
    #if we are running the fourth experiment    
    elif Etype == 2:
        
        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~') 
        cwd = os.getcwd()
        #run the ADP algorithm
        ADP_dict, ADP_alleles, ADP_prop = run_ADP(editDist, iteration)
        strain_df = pf.strainSolver(cwd,strainRef,outputDir,genes,'noProp','s',10800,5)
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        ADP_alleles = set(item for sublist in ADP_alleles for item in sublist)
        #compute the precision and recall for the SDP output
        ADP_pred, ADP_recall, ADP_tvd = Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        return prec, rec, tvd, ADP_pred, ADP_recall, ADP_tvd
        
    
def Process_reads(gene, editDist, iteration):
    mapping_cmd = "bowtie -a --best --strata -v 3 -p 4 {0}_bowtie -1 ./editDist_{1}_iteration_{2}_all_Strains_1.fa -2 ./editDist_{1}_iteration_{2}_all_Strains_2.fa --sam all_Strains.sam >/dev/null 2>&1".format(gene, editDist, iteration)
    os.system(mapping_cmd)
    mapped_cmd = "samtools view -h -F4 all_Strains.sam > all_Strains_mapped.sam"
    paired_cmd = "samtools view -F8 all_Strains_mapped.sam > all_Strains_pairedNoHeader.sam"
    os.system(mapped_cmd)
    os.system(paired_cmd)
    writeReadTable()
    pairedDF = pf.returnMismatchMatrix('all_Strains_reads.txt', "paired")
    dataMatrix = pairedDF
    dataMatrix = dataMatrix.fillna(-1)
    paired_Qmatrix = pf.returnQualityMatrix('all_Strains_reads.txt', "paired")
#       singleton_Qmatrix = returnQualityMatrix(singleton_readsTxt_path, "singleton")
    Qmatrix = paired_Qmatrix
    
        #Run the ILP solver
    pred_object_val,var_predicted,reads_cov,all_solutions, all_objective = varSolver.solver(dataMatrix)
    predicted_DF = dataMatrix.loc[reads_cov,var_predicted]
    prop_count = pf.compute_proportions(predicted_DF)
        #prop_bayes = bayes_compute_proportions(predicted_DF)
#        if proportion_method == "count":
#            prop = count_compute_proportions(predicted_DF)
#        elif proportion_method == "bayes":
#            prop = bayes_compute_proportions(predicted_DF)
    pred_prop_count = pf.create_dictionary(var_predicted, prop_count)
    return pred_prop_count

#A function to run the ADP on individual genes
def run_ADP(editDist, iteration):
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    ADP_dict = defaultdict(int)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE ALLELE DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    ADP_alleles = []
    ADP_prop = []
    for gene in genes:
        print('now processing gene: {}'.format(gene))
        result = Process_reads(gene, editDist, iteration)
        for key in result.keys():
            if key not in ADP_dict:
                ADP_dict[key] = result[key]
            else:
                ADP_dict[key] = ADP_dict[key] + result[key]
        #os.chdir(samplesDir)
        Dict_to_csv(result,gene)
        ADP_alleles.append(result.keys())
        ADP_prop.append(result.values())
        #print('The predicted genes and their proportions are: {}'.format(result))
        #print('The sum of the alleles proportions are: {}'.format(sum(result.values())))
        #os.chdir('../')
        #ADP_dict[gene] = result
    return ADP_dict, ADP_alleles, ADP_prop
    
def main():
    
    #get and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--numOfiter", required=False, default=40, type=int, help="The number of simulation iterations default = 40")
    ap.add_argument("-d", "--masterDir", required=False, default='Results',type=str, help="The name of the master directory to run the simulation and store the results. Default = Created folder called Results")
    #ap.add_argument("-r", "--strainRef", required=True, help="Name of the text file containing the reference strains that is stored in this directory")
    #ap.add_argument("-t", "--Sim_Type", required=False, default='Mutation', type=str, help="The type of simulation you would like to run (Mutation or Recombination). Default = Mutation")
    ap.add_argument("-Ed", "--editDist", required=False, default=5, type=int, help="The maximum edit distance. Default = 5")
    ap.add_argument("-nm", "--numMut", required=False, default=1, type=int, help="The number of mutations to introduce. Default = 1")
    ap.add_argument("-st", "--modType", required=False, default="Type_1", type=str, help="The type of evolution experiment to run. Default = Type 1")
    ap.add_argument("-et", "--experimentType", required=False, default=1, type=int, help="The type of experiment to run. Default = 1")
    args = vars(ap.parse_args())
    
    directories = [d for d in os.listdir(".") if os.path.isdir(d)]
    if args['masterDir'] not in directories:
        os.mkdir(args['masterDir'])
    
    precision = []
    recall = []
    total_var_dist = []
    ADP_prec_vals = []
    ADP_rec_vals = []
    ADP_tvd_vals = []
    if args["experimentType"] == 1:
        if args["modType"] == "Type_1":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd = EvoMod1(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
        elif args["modType"] == "Type_2":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd = EvoMod2(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
        elif args["modType"] == "Type_3":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd = EvoMod3(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
             
    elif args["experimentType"] == 2:
        if args["modType"] == "Type_1":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, ADP_prec, ADP_rec, ADP_tvd = EvoMod1(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                ADP_prec_vals.append(ADP_prec)
                ADP_rec_vals.append(ADP_rec)
                ADP_tvd_vals.append(ADP_tvd)
        elif args["modType"] == "Type_2":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, ADP_prec, ADP_rec, ADP_tvd = EvoMod2(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                ADP_prec_vals.append(ADP_prec)
                ADP_rec_vals.append(ADP_rec)
                ADP_tvd_vals.append(ADP_tvd)
        elif args["modType"] == "Type_3":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, ADP_prec, ADP_rec, ADP_tvd = EvoMod3(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                ADP_prec_vals.append(ADP_prec)
                ADP_rec_vals.append(ADP_rec)
                ADP_tvd_vals.append(ADP_tvd)
            
    
    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    avg_prec = sum(precision)/len(precision)
    std_prec = np.std(np.array(precision))
    print 'Precision is:', precision, "\n"
    print "Average of precision is: ", avg_prec, "\n"
    print "Standard deviation of precision is: ", std_prec, "\n"

    
    avg_rec = sum(recall)/len(recall)
    std_rec = np.std(np.array(recall))
    print 'Recall is:', recall, "\n"
    print "Average of recall is: ", avg_rec, "\n"
    print "Standard deviation of recall is: ", std_rec, "\n"
    
    avg_TVD = sum(total_var_dist)/len(total_var_dist)
    std_TVD = np.std(np.array(total_var_dist))
    print 'Total Variation Distance is:', total_var_dist, "\n"
    print "Average of Total Variation Distance is: ", avg_TVD, "\n"
    print "Standard deviation of Total Variation Distance is: ", std_TVD, "\n"
    
    if args["experimentType"] == 2:
        os.chdir(args["masterDir"])
        
        avg_ADP_prec = sum(ADP_prec_vals)/len(ADP_prec_vals)
        std_ADP_prec = np.std(np.array(ADP_prec_vals))
        print 'Precision for ADP is:', ADP_prec_vals, "\n"
        print "Average of ADP precision is: ", avg_ADP_prec, "\n"
        print "Standard deviation of ADP precision is: ", std_ADP_prec, "\n"
        
        ADPprecisionDF = pd.DataFrame(ADP_prec_vals)
        ADPprecisionDF = ADPprecisionDF.T
        ADPprecisionDF.to_csv('ADP_{}_editDist_{}_Precision_values.csv'.format(args["modType"], args["editDist"]),sep = '\t')

    
        avg_ADP_rec = sum(ADP_rec_vals)/len(ADP_rec_vals)
        std_ADP_rec = np.std(np.array(ADP_rec_vals))
        print 'Recall for ADP is:', ADP_rec_vals, "\n"
        print "Average of ADP recall is: ", avg_ADP_rec, "\n"
        print "Standard deviation of ADP recall is: ", std_ADP_prec, "\n"
        
        ADPrecallDF = pd.DataFrame(ADP_rec_vals)
        ADPrecallDF = ADPrecallDF.T
        ADPrecallDF.to_csv('ADP_{}_editDist_{}_Recall_values.csv'.format(args["modType"],args["editDist"]), sep = '\t')
    
        ADP_avg_TVD = sum(ADP_tvd_vals)/len(ADP_tvd_vals)
        ADP_std_TVD = np.std(np.array(ADP_tvd_vals))
        print 'ADP Total Variation Distance is:', ADP_tvd_vals, "\n"
        print "Average of ADP Total Variation Distance is: ", ADP_avg_TVD, "\n"
        print "Standard deviation of ADP Total Variation Distance is: ", ADP_std_TVD, "\n"
        
        ADPtotal_var_distDF = pd.DataFrame(ADP_tvd_vals)
        ADPtotal_var_distDF = ADPtotal_var_distDF.T
        ADPtotal_var_distDF.to_csv('ADP_{}_editDist_{}_Total_Variation_Distance.csv'.format(args["modType"], args["editDist"]), sep = '\t')
        
        os.chdir('../')
    
    
    os.chdir(args["masterDir"])
    plt.figure()
    plt.hist(recall, bins=np.linspace(0,2))
    plt.title("Plot of simulation recall for {} for edit Distance {}".format(args["modType"], args["editDist"]))
    plt.xlabel("Recall")
    plt.ylabel("Frequency")
    plt.savefig('{}_editDist_{}_Recall_plot'.format(args["modType"],args["editDist"]))
    #Save the plot for boxplot plotting
    recallDF = pd.DataFrame(recall)
    recallDF = recallDF.T
    recallDF.to_csv('{}_editDist_{}_Recall_values.csv'.format(args["modType"],args["editDist"]), sep = '\t')
    #for precision
    plt.figure()
    plt.hist(precision, bins=np.linspace(0,2))
    plt.title("Plot of simulation precision for {} for edit Distance {}".format(args["modType"], args["editDist"]))
    plt.xlabel("Precision")
    plt.ylabel("Frequency")
    plt.savefig('{}_editDist_{}_Precision_plot'.format(args["modType"], args["editDist"]))
    #save the plot for boxplot plotting
    precisionDF = pd.DataFrame(precision)
    precisionDF = precisionDF.T
    precisionDF.to_csv('{}_editDist_{}_Precision_values.csv'.format(args["modType"], args["editDist"]),sep = '\t')
    #for total variation distance
    plt.figure()
    plt.hist(total_var_dist,bins=np.linspace(-1,1))
    plt.title("Plot of simulation Total variation distance for {} for edit Distance {}".format(args["modType"], args["editDist"]))
    plt.xlabel("Total Variation Distance")
    plt.ylabel("Frequency")
    plt.savefig('{}_editDist_{}_Total_variation_Distance_plot'.format(args["modType"], args["editDist"]))
    #save the plot for boxplot plotting
    total_var_distDF = pd.DataFrame(total_var_dist)
    total_var_distDF = total_var_distDF.T
    total_var_distDF.to_csv('{}_editDist_{}_Total_Variation_Distance.csv'.format(args["modType"], args["editDist"]), sep = '\t')
    
    
#run the main function
if __name__ == "__main__":
    main()
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    