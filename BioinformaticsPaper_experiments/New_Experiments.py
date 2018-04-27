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
import sys
import os
plt.style.use('ggplot')

#Global variable
EDITDIST = 5    #For soft precision and recall calculation
MAX_WEIGHT = 100

#Set seed for reproducibility
seed = 1995
random.seed(seed)

#set up genes edit distance matrices

project_path = os.getcwd()

#clpA
clpA_df = pd.read_csv(os.path.join(project_path,'editDistanceMatrix_clpA.csv'), sep=",")
#clpX
clpX_df = pd.read_csv(os.path.join(project_path, 'editDistanceMatrix_clpX.csv'), sep=",")
#nifS
nifS_df = pd.read_csv(os.path.join(project_path, 'editDistanceMatrix_nifS.csv'), sep=",")
#pepX
pepX_df = pd.read_csv(os.path.join(project_path, 'editDistanceMatrix_pepX.csv'), sep=",")
#pyrG
pyrG_df = pd.read_csv(os.path.join(project_path, 'editDistanceMatrix_pyrG.csv'), sep=",")
#recG
recG_df = pd.read_csv(os.path.join(project_path, 'editDistanceMatrix_recG.csv'), sep=",")
#rplB
rplB_df = pd.read_csv(os.path.join(project_path, 'editDistanceMatrix_rplB.csv'), sep=",")
#uvrA
uvrA_df = pd.read_csv(os.path.join(project_path, 'editDistanceMatrix_uvrA.csv'), sep=",")

#read in the reference strain database
loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
genes_df = [clpA_df,clpX_df,nifS_df,pepX_df,pyrG_df,recG_df,rplB_df,uvrA_df]
genes_dict = {locus:[genes_df[i]] for i,locus in enumerate(loci)}
reference = pd.read_csv(os.path.join(project_path, 'New_strain_ref.txt'),sep="\t",usecols=range(1,len(loci)+1))
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


def compute_likelihood(df, max_mm):
    numVar = df.shape[1]
    likelihood_list = list()

    for row in df.itertuples(index=False):
        read = list(row)

        temp = list()
        for i in range(numVar):
            if read[i] == -1:   #treat those reads which do not map having mm=max_mm+1
                prob = (0.01)**(max_mm+1) * (0.99)**(152 - max_mm -1)
                temp.append(prob)
            else:
                prob = (0.01)**(read[i]) * (0.99)**(152 - read[i])
                temp.append(prob)

        likelihood_list.append( sum(temp) )

    #Similar to method in GAML paper
    likelihood_list = [i/(2.0*152*numVar) for i in likelihood_list]
    neg_log_likelihood = [-1.0*np.log10(j) for j in likelihood_list]

    score = sum(neg_log_likelihood)
    return score

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


def Generate_reads(strains, seed, iteration, editDist, art_cmd):
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
        ART_command = art_cmd + " -qL 33 -qs 10 -qs2 15 -k 3 -rs {0} -q -ss HS25 -sam -i {1} -p -l 76 -f {2} -m 149.45 -s 41.35 -o {3}".format(seed, 'editDist_{}_iteration_{}_Strain_{}_sequences.fas >/dev/null 2>&1'.format(editDist,iteration,i),proportions[i]*100, output_file_name)
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

def Compute_Soft_Prec_and_Rec(strain_dict, strain_df, editDist):
    '''
        A function that returns soft precision and recall
        strain_dict: A dictionary containing information about true strains where
                    keys(a list, the strains), values are proportions
        strain_df: A dataframe where there are ST, New/Existing and 8 loci as columns, and Proportions column
    '''

    #False positive
    #FP = Predicted set of trains - True set of strains
    #i/j is a list of alleles, sort it so that we are comparing correctly after converting to tuples
    #Use set operations to get the false positives
    PredictedStrains = [ tuple(sorted(i)) for i in strain_df.values[:,2:-1] ]
    TrueStrains = [ tuple(sorted(i)) for i in strain_dict.keys() ]
    FP = list( set(PredictedStrains) - set( TrueStrains ) )
    #True positives
    TP = [i for i in PredictedStrains if i in TrueStrains]

    #False negatives
    FN = [i for i in TrueStrains if i not in PredictedStrains]

    #Number of soft true positives i.e. number of real true positives + strains which are in acceptable range
    numSoft_TP = len(TP)
    for fp in FP:
        fpIsMatched = False     #indicator if a false positive is matched

        trueS_track = 0
        while( (trueS_track <= len(TrueStrains)-1) and (fpIsMatched == False) ):
            sumOfEd = 0     #The edit distance for a whole strain

            trueS = TrueStrains[trueS_track]

            #Iterate over each locus
            for locus in range(len(trueS)):
                #If they are the same, pass
                if fp[locus] == trueS[locus]:
                    pass
                else:
                    loc_name = trueS[locus].split('_')[0]

                    #hard coded
                    #Retrieve edit distance
                    if loc_name == 'clpA':
                        ed = int(clpA_df.loc[clpA_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed
                    elif loc_name == 'clpX':
                        ed = int(clpX_df.loc[clpX_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed
                    elif loc_name == 'nifS':
                        ed = int(nifS_df.loc[nifS_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed
                    elif loc_name == 'pepX':
                        ed = int(pepX_df.loc[pepX_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed
                    elif loc_name == 'pyrG':
                        ed = int(pyrG_df.loc[pyrG_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed
                    elif loc_name == 'recG':
                        ed = int(recG_df.loc[recG_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed
                    elif loc_name == 'rplB':
                        ed = int(rplB_df.loc[rplB_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed
                    elif loc_name == 'uvrA':
                        ed = int(uvrA_df.loc[uvrA_df['level_0'] == fp[locus]][trueS[locus]])
                        sumOfEd += ed

            #If edit distance as a strain is in acceptable range
            #Accept it as true
            if sumOfEd <= editDist:
                numSoft_TP += 1
                fpIsMatched = True

            trueS_track += 1
    soft_precision = float(numSoft_TP)/float(len(PredictedStrains))
    soft_recall = float(numSoft_TP)/len(TrueStrains)

    #Makes sure it's max is 1.0
    if soft_precision > 1.0:
        soft_precision = 1.0
    if soft_recall > 1.0:
        soft_recall = 1.0


    return soft_precision, soft_recall

'''
Returns the edit distance between two alleles
Input:
        allele_p: predicted allele, a string
        allele_t: true_allele, a string
        gene: a string
'''

def compute_distance_by_gene(allele_p, allele_t,gene):
    assert(gene in genes_dict.keys()), "{} not in the list of genes".format(gene)
    assert(allele_p in genes_dict[gene][0].columns), "{} not in the allele list".format(allele_p)
    assert(allele_t in genes_dict[gene][0].columns), "{} not in the allele list".format(allele_t)
    return int( genes_dict[gene][0].loc[ genes_dict[gene][0]['level_0'] == allele_p ][allele_t] )

'''
Computes the distance between two strains. It is assumed that pred and true has the same length and each position represents the same gene.
Input:
        pred: a list, containing an allele at each gene/position. Represents predicted strains
        true: '' Represents true strains
        mode: 'snp' computes edit distances in terms of bases, 'hamming' computes hamming distance in terms of gene
'''

def compute_strain_distance(pred, true, mode="snp"):
    assert(len(pred) == len(true)), "pred and true needs to be the same length"
    assert(mode == "snp" or mode == "hamming"), "Only 'snp' or 'hamming' options are permitted"

    if mode == "hamming":
        from scipy.spatial.distance import hamming
        dist = hamming(pred, true)
    else:
        strain_length = len(pred)
        
        dist = 0
        for i in range(strain_length):
            p_gene = pred[i].split("_")[0]
            t_gene = true[i].split("_")[0]

            assert(p_gene == t_gene), "pred and true have different genes at position {}".format(i)
            dist += compute_distance_by_gene(pred[i], true[i], p_gene)

    return dist


'''
Returns the Earth-Mover's distance between predicted strains and true strains.
Input:
        pred_dict: A dictionary for predicted strains, where key=strain(tuple/list of multi-loci), value=proportion of the strain
        true_dict: A dictionary for true strains, where ...
        mode: 'hamming' or 'snp'. 'hamming' computes the distance in terms of gene and 'snp' is in terms bases.
'''

def compute_EMD(pred_dict, true_dict, mode="snp"):
    #Check if pred_list and true_list are dictionaries
    #assert(type(pred_dict) is dict), "pred_dict is not a dictionary"
    #assert(type(true_dict) is dict), "true_dictionary is not a dictionary"

    #Convert to lists.
    pred_list = pred_dict.keys()
    true_list = true_dict.keys()

    #Makes sure gene are comparable. Convert to list if it isn't
    #pred_list = [sorted(list(i)) for i in pred_list]
    #true_list = [sorted(list(i)) for i in true_list]

    #Create a dictionary for keeping track strains
    pred_s_dict = {"p{}".format(i): strain for i,strain in enumerate(pred_list)}
    true_s_dict = {"t{}".format(i): strain for i,strain in enumerate(true_list)}

    #Short form 
    pred_s_list = np.array(pred_s_dict.keys())
    true_s_list = np.array(true_s_dict.keys())

    #Convert proportions to numpy array, it serves as an input to emd weights
    pred_prop = np.array([pred_dict[ tuple(pred_s_dict[s]) ] for s in pred_s_list])
    true_prop = np.array([true_dict[ tuple(true_s_dict[s]) ] for s in true_s_list])

    distance_dict = {}
    distance_mat = np.zeros((pred_s_list.shape[0], true_s_list.shape[0]))

    for i,p in enumerate(pred_s_list):
        for j,t in enumerate(true_s_list):
            if pred_s_dict[p] == true_s_dict[t]:
                distance_dict[(p,t)] = 0
                distance_mat[i][j] = 0
            else:
                d = compute_strain_distance(pred_s_dict[p],true_s_dict[t],mode)
                distance_dict[(p,t)] = d
                distance_mat[i][j] = d

    #Compute emd
    from emd import emd
    #emd function only accepts np array
    pred_s_array = np.array(pred_s_list).reshape(len(pred_s_list),-1)
    true_s_array = np.array(true_s_list).reshape(len(true_s_list), -1)

    EMD = emd(pred_s_array, true_s_array,X_weights=pred_prop, Y_weights=true_prop, distance='precomputed', D=distance_mat)
    #import pdb
    #pdb.set_trace()
    return EMD 
    
    

def Compute_EditDist_Stats(strain_dict, strain_df):
    '''
    This function tries to map each predicted strain to the closest 1 true strain.
    The weight of the edge of this mapping is the edit distance between the strains, which the sum of
    the weight of all edges are then returned.
    If a true strain is not paired with any predicted strain, we penalize a high penalty defined by MAX_WEIGHT.

    Returns: weight, sum of weightage of the edges for bipartite matching
            numRand, number of times when the randomized section of the algorithm is utilized
    strain_dict: A dictionary containing information about true strains where
                keys(a list, the strains), values are proportions
    strain_df: A dataframe where there are ST, New/Existing and 8 loci as columns, and Proportions column
'''
    #False positive
    #FP = Predicted set of trains - True set of strains
    #i/j is a list of alleles, sort it so that we are comparing correctly after converting to tuples
    #Use set operations to get the false positives
    PredictedStrains = [ tuple(sorted(i)) for i in strain_df.values[:,2:-1] ]
    TrueStrains = [ tuple(sorted(i)) for i in strain_dict.keys() ]
    FP = list( set(PredictedStrains) - set( TrueStrains ) )
    #True positives
    TP = [i for i in PredictedStrains if i in TrueStrains]

    #False negatives
    FN = [i for i in TrueStrains if i not in PredictedStrains]

    #Create random floats between 0 and 1 for randomized algorithm
    randTrue = np.random.uniform(0,1,len(TrueStrains))
    randTrue_dict = {TrueStrains[i]: randTrue[i] for i in range(len(TrueStrains))}

    #sum of weights to report
    weight = 0

    #Count the number of matchings for each true strains
    #For true positives, the match is automatically set to 1 for now
    count_TrueS = {i:0 for i in TrueStrains}

    #For true positives, the match is automatically set to 1 for now
    for i in TP:
        count_TrueS[i] = 1

    #Keep track the mapping of FP to which true strain
    track_FP = {i: None for i in FP}

    #Number of times random section of the algo is utilized
    numRand = 0

    #A very greedy randomized algorithm for matching
    for fp in FP:
        ed_list = list()

        for ts in TrueStrains:
            sumOfEd = 0
            for locus in range(len(ts)):
                #If same allele
                if fp[locus] == ts[locus]:
                    pass
                else:
                    loc_name = ts[locus].split('_')[0]

                    #hard coded
                    #Retrieve edit distance
                    if loc_name == 'clpA':
                        ed = int(clpA_df.loc[clpA_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed
                    elif loc_name == 'clpX':
                        ed = int(clpX_df.loc[clpX_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed
                    elif loc_name == 'nifS':
                        ed = int(nifS_df.loc[nifS_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed
                    elif loc_name == 'pepX':
                        ed = int(pepX_df.loc[pepX_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed
                    elif loc_name == 'pyrG':
                        ed = int(pyrG_df.loc[pyrG_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed
                    elif loc_name == 'recG':
                        ed = int(recG_df.loc[recG_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed
                    elif loc_name == 'rplB':
                        ed = int(rplB_df.loc[rplB_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed
                    elif loc_name == 'uvrA':
                        ed = int(uvrA_df.loc[uvrA_df['level_0'] == fp[locus]][ts[locus]])
                        sumOfEd += ed

            #Keep track of all edit distances with true strains
            ed_list.append(sumOfEd)

        #Find the minimum edit distance
        min_ed = np.min(ed_list)

        #Find the strains which has this minimum
        min_ed_trueS_index = [i for i in range(len(ed_list)) if ed_list[i] == min_ed ]
        min_ed_trueS = [TrueStrains[s] for s in min_ed_trueS_index]

        #Increment the weight by minimum edit dist, doesn't depend on the following algo
        weight += min_ed

        #If there is only one minimum strain
        if len(min_ed_trueS) == 1:
            count_TrueS[min_ed_trueS[0]] += 1
            track_FP[fp] = min_ed_trueS[0]

        #If there are more than one minimum strains
        else:
        #match to the one with smallest degree
            min_deg = float('inf')

            for s in min_ed_trueS:
                if count_TrueS[s] < min_deg:
                    min_deg = count_TrueS[s]

            min_deg_str = [s for s in min_ed_trueS if count_TrueS[s] == min_deg]

        #If there is only one min degree strain, match to that
            if len(min_deg_str) == 1:
                count_TrueS[min_deg_str[0]] += 1
                track_FP[fp] = min_deg_str[0]
            else:
        #randomized algorithm for matching, assign to minimum float
                numRand += 1
                min_rand = float("inf")
                finalMatch = None

                for s in min_deg_str:
                    if randTrue_dict[s] < min_rand:
                        min_rand = randTrue_dict[s]
                        finalMatch = s

                count_TrueS[finalMatch] += 1
                track_FP[fp] = finalMatch

    #Check if any true strains are not mapped yet. Suffices to check those in FN
    unmapped_FN = [i for i in FN if count_TrueS[i] == 0]

    #Penalize weight with high penality defined by MAX_WEIGHT if there are FN unmapped
    weight += MAX_WEIGHT*len(unmapped_FN)

    print("\n*****************")
    print("True strains are: {}\n".format(TrueStrains))
    print("Predicted strains are: {}\n".format(PredictedStrains))
    print("False positive mappings: {}\n".format(track_FP))
    print("Number of randomized event occurred: {}\n".format(numRand))
    print("Degree of nodes representing true strains: {}\n".format(count_TrueS))
    print("*****************\n")

    return weight, numRand

def Compute_Prec_and_rec(strain_dict, strain_df):
    '''
        A function to compute the precision and recall after a simulation.

        strain_dict: A dictionary containing the true strains

        strain_df: A dataframe containing the strains predicted by the ILP

    '''
    #import pdb
    #pdb.set_trace()
    count = 0
    total_var_dist = 0
    predicted_dict = defaultdict(list)
    #convert the dataframe into a dictionary
    for i in range(strain_df.shape[0]):
        row = strain_df.iloc[i].tolist()
        temp_key = tuple(row[2:10])
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
    #pdb.set_trace()

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

#A function to compute the set of known strains in the library without the
#new strains to the set of predicted strains.
def CompareKnown(true, predicted):
    count = 0
    total_var_dist = 0
    predicted_dict = defaultdict(list)
    #convert the dataframe into a dictionary
    for i in range(predicted.shape[0]):
        row = predicted.iloc[i].tolist()
        if row[1] == 'Existing':
            temp_key = tuple(row[2:10])
            predicted_dict[temp_key] = float(row[10])
    print ('The Known strains and their proportions are {}\n'.format(true))
    print ('The Predicted strains and their proportions are\n {}'.format(predicted_dict))
    #Compute the total variation distance
    for key in predicted_dict.keys():
        if key in true.keys():
            total_var_dist += abs(round(true[key],3)-predicted_dict[key])
            count += 1
        else:
            total_var_dist += predicted_dict[key]
    for key in true.keys():
	if key not in predicted_dict.keys():
	    total_var_dist += true[key]
    #compute the precision
    try:
        precision = count/len(predicted_dict)
    except ZeroDivisionError:
        precision = 0

    #compute the recall
    try:
        recall = count/len(true)
    except ZeroDivisionError:
        recall = 0
        
    return precision, recall, total_var_dist/2.0



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
def EvoMod1(reference, editDist, num_mut, iteration, Type, Etype, distance_mode='snp',art_cmd="art_illumina"):
    project_path = os.getcwd()
    strainRef = os.path.join(project_path,'New_strain_ref.txt')
    samplesDir = os.path.join(project_path, 'Type_1_editDist_{}_iteration_{}_SDP_Samples'.format(editDist,iteration))
    outputDir = os.path.join(project_path, 'Type_1_editDist_{}_iteration_{}_SDP_output'.format(editDist,iteration))
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
    proportions = Generate_reads(strains, seed, iteration, editDist, art_cmd)
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
        true_dict = strain_prop_dict
        pred_dict = dict()
        for i in range(strain_df.shape[0]):
            row = strain_df.iloc[i].tolist()
            temp_key = tuple(row[2:10])
            pred_dict[temp_key] = float(row[10])

        #Sort the genes in the strains to be consistent
        true_dict = {tuple(sorted(i)):true_dict[i] for i in true_dict.keys()}
        pred_dict = {tuple(sorted(i)):pred_dict[i] for i in pred_dict.keys()}
        #compute the precision and recall
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        soft_prec, soft_rec = Compute_Soft_Prec_and_Rec(strain_prop_dict, strain_df, EDITDIST)
        EMD = compute_EMD(pred_dict, true_dict, distance_mode)
        #weight, numRand = Compute_EditDist_Stats(strain_prop_dict, strain_df)
        return prec, rec, tvd, soft_prec, soft_rec, EMD
    #if we are running the fourth experiment
    elif Etype == 2:

        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        cwd = os.getcwd()
        #run the ADP algorithm
        ADP_dict, ADP_alleles, ADP_prop = run_ADP(editDist, iteration)
        strain_df = pf.strainSolver(cwd,strainRef,outputDir,genes,'noProp','s',10800,5)
        true_dict = strain_prop_dict
        pred_dict = dict()
        for i in range(strain_df.shape[0]):
            row = strain_df.iloc[i].tolist()
            temp_key = tuple(row[2:10])
            pred_dict[temp_key] = float(row[10])
        #Sort the genes in the strains to be consistent
        true_dict = {tuple(sorted(i)):true_dict[i] for i in true_dict.keys()}
        pred_dict = {tuple(sorted(i)):pred_dict[i] for i in pred_dict.keys()}
        SDP_list = []
        ADP_alleles = set(item for sublist in ADP_alleles for item in sublist)
        #compute the precision and recall for the SDP output
        ADP_pred, ADP_recall, ADP_tvd = Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        soft_prec, soft_rec = Compute_Soft_Prec_and_Rec(strain_prop_dict, strain_df, EDITDIST)
        EMD = compute_EMD(pred_dict, true_dict, distance_mode)
        #weight, numRand = Compute_EditDist_Stats(strain_prop_dict, strain_df)
        return prec, rec, tvd, soft_prec, soft_rec, EMD, ADP_pred, ADP_recall, ADP_tvd

'''
EvoMod2 (Result in 4 strains, 2 new 2 existing)
1) Pick 2 strains at random from database
2) Repeat EvoMod1 on each of these strains, make sure that the 2 mutated strains are novel strains
3) Now we have 4 strains
4) Same stuff as above
'''

def EvoMod2(reference, editDist, num_mut, iteration, Type, Etype, distance_mode='snp',art_cmd="art_illumina"):
    project_path = os.getcwd()
    strainRef = os.path.join(project_path,'New_strain_ref.txt')
    samplesDir = os.path.join(project_path, 'Type_2_editDist_{}_iteration_{}_SDP_Samples'.format(editDist,iteration))
    outputDir = os.path.join(project_path, 'Type_2_editDist_{}_iteration_{}_SDP_output'.format(editDist,iteration))
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    strains = []
     #get the indices of the strains to use
    indices = random.sample(range(1, 696), 2)
    known_strains = defaultdict(int)
    for index in indices:
        S1 = reference.iloc[index,:].values.tolist()
        strains.append(S1)
        #add the known strains to dictionary
        known_strains[tuple(S1)] = 1
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
    proportions = Generate_reads(strains, seed, iteration, editDist, art_cmd)
    for strain in strains:
        #if this is a strain already in the database, update its proportions
        if strain in known_strains.keys():
            known_strains[strain] = proportions[strains.index(strain)]
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
        true_dict = strain_prop_dict
        pred_dict = dict()
        for i in range(strain_df.shape[0]):
            row = strain_df.iloc[i].tolist()
            temp_key = tuple(row[2:10])
            pred_dict[temp_key] = float(row[10])
        #Sort the genes in the strains to be consistent
        true_dict = {tuple(sorted(i)):true_dict[i] for i in true_dict.keys()}
        pred_dict = {tuple(sorted(i)):pred_dict[i] for i in pred_dict.keys()}
        #compute the precision and recall
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        soft_prec, soft_rec = Compute_Soft_Prec_and_Rec(strain_prop_dict, strain_df, EDITDIST)
        EMD = compute_EMD(pred_dict, true_dict, distance_mode)
        #weight, numRand = Compute_EditDist_Stats(strain_prop_dict, strain_df)
        return prec, rec, tvd, soft_prec, soft_rec, EMD
    #if we are running the fourth experiment
    elif Etype == 2:

        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        cwd = os.getcwd()
        #run the ADP algorithm
        ADP_dict, ADP_alleles, ADP_prop = run_ADP(editDist, iteration)
        strain_df = pf.strainSolver(cwd,strainRef,outputDir,genes,'noProp','s',10800,5)
        true_dict = strain_prop_dict
        pred_dict = dict()
        for i in range(strain_df.shape[0]):
            row = strain_df.iloc[i].tolist()
            temp_key = tuple(row[2:10])
            pred_dict[temp_key] = float(row[10])
        #Sort the genes in the strains to be consistent
        true_dict = {tuple(sorted(i)):true_dict[i] for i in true_dict.keys()}
        pred_dict = {tuple(sorted(i)):pred_dict[i] for i in pred_dict.keys()}
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        soft_prec, soft_rec = Compute_Soft_Prec_and_Rec(strain_prop_dict, strain_df, EDITDIST)
        ADP_alleles = set(item for sublist in ADP_alleles for item in sublist)
        #compute the precision and recall for the SDP output
        ADP_pred, ADP_recall, ADP_tvd = Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        #weight, numRand = Compute_EditDist_Stats(strain_prop_dict, strain_df)
        EMD = compute_EMD(pred_dict, true_dict, distance_mode)
        known_prec, known_rec, known_tvd = CompareKnown(known_strains, strain_df)
        return prec, rec, tvd, soft_prec, soft_rec, EMD, ADP_pred, ADP_recall, ADP_tvd, known_prec, known_rec


'''
EvoMod3(Result in 5 strains, 3 new 2 existing)
1) Do EvoMod2
2) Choose any 2 strains to recombine, to produce a 3rd new strain
3) If we can afford to do, 2 recombination events. We could try that
4) ...
'''

def EvoMod3(reference, editDist, num_mut, iteration, Type, Etype, distance_mode='snp', art_cmd="art_illumina"):
    #do EveMod2
    project_path = os.getcwd()
    strainRef = os.path.join(project_path,'New_strain_ref.txt')
    samplesDir = os.path.join(project_path, 'Type_3_editDist_{}_iteration_{}_SDP_Samples'.format(editDist,iteration))
    outputDir = os.path.join(project_path, 'Type_3_editDist_{}_iteration_{}_SDP_output'.format(editDist,iteration))
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
    proportions = Generate_reads(strains, seed, iteration, editDist, art_cmd)
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
        true_dict = strain_prop_dict
        pred_dict = dict()
        for i in range(strain_df.shape[0]):
            row = strain_df.iloc[i].tolist()
            temp_key = tuple(row[2:10])
            pred_dict[temp_key] = float(row[10])
        #Sort the genes in the strains to be consistent
        true_dict = {tuple(sorted(i)):true_dict[i] for i in true_dict.keys()}
        pred_dict = {tuple(sorted(i)):pred_dict[i] for i in pred_dict.keys()}
        #compute the precision and recall
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        soft_prec, soft_rec = Compute_Soft_Prec_and_Rec(strain_prop_dict, strain_df, EDITDIST)
        EMD = compute_EMD(pred_dict, true_dict, distance_mode)
        #weight, numRand = Compute_EditDist_Stats(strain_prop_dict, strain_df)
        return prec, rec, tvd, soft_prec, soft_rec, EMD
    #if we are running the fourth experiment
    elif Etype == 2:

        strain_prop_dict = Write_Proportions(strains, proportions, iteration, editDist, Type)
        print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAIN DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        cwd = os.getcwd()
        #run the ADP algorithm
        ADP_dict, ADP_alleles, ADP_prop = run_ADP(editDist, iteration)
        strain_df = pf.strainSolver(cwd,strainRef,outputDir,genes,'noProp','s',10800,5)
        true_dict = strain_prop_dict
        pred_dict = dict()
        for i in range(strain_df.shape[0]):
            row = strain_df.iloc[i].tolist()
            temp_key = tuple(row[2:10])
            pred_dict[temp_key] = float(row[10])
        #Sort the genes in the strains to be consistent
        true_dict = {tuple(sorted(i)):true_dict[i] for i in true_dict.keys()}
        pred_dict = {tuple(sorted(i)):pred_dict[i] for i in pred_dict.keys()}
        SDP_list = []
        for row in strain_df.iterrows():
            index, data = row
            SDP_list.append(data.tolist()[0:8])
        prec, rec, tvd = Compute_Prec_and_rec(strain_prop_dict, strain_df)
        soft_prec, soft_rec = Compute_Soft_Prec_and_Rec(strain_prop_dict, strain_df, EDITDIST)
        #weight, numRand = Compute_EditDist_Stats(strain_prop_dict, strain_df)
        ADP_alleles = set(item for sublist in ADP_alleles for item in sublist)
        #compute the precision and recall for the SDP output
        ADP_pred, ADP_recall, ADP_tvd = Compute_ADP_Prec_and_rec(true_all, ADP_dict)
        return prec, rec, tvd, soft_prec, soft_rec, EMD, ADP_pred, ADP_recall, ADP_tvd

def upperFirst(x):
    return x[0].upper() + x[1:]

def Process_reads(gene, editDist, iteration):
    mapping_cmd = "bowtie -a --best --strata -v 3 -p 4 {0}_bowtie -1 ./editDist_{1}_iteration_{2}_all_Strains_1.fa -2 ./editDist_{1}_iteration_{2}_all_Strains_2.fa --sam all_Strains.sam".format(upperFirst(gene), editDist, iteration)
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
    pred_object_val,var_predicted,reads_cov,all_solutions, all_objective = varSolver.solver(dataMatrix, paired_Qmatrix, None)
    #compute the likelihood scores
    score_list = list()
    min_score = sys.maxint

    #Compute negative log likelihood score for each solution
    for i in range(len(all_solutions)):
        score = pf.compute_likelihood(dataMatrix.loc[reads_cov, all_solutions[i]], 6)
        score_list.append(score)

        if score <= min_score:
            min_score = score

    argmin_score_list = [i for i in range(len(all_solutions)) if score_list[i] == min_score]
    if len(argmin_score_list) > 1:
        print("More than 1 solution having minimum negative log likelihood score.")
        lexico_min_score_sol = [all_solutions[i] for i in argmin_score_list]
        lexico_min_score_sol = sorted(lexico_min_score_sol)
        var_predicted = lexico_min_score_sol[0]
    else:
        var_predicted = all_solutions[argmin_score_list[0]]
    print("The minimum negative log likelihood score(s) are: {}".format(score_list))
    print("The chosen solution based on minimum negative log likelihood is: {}".format(var_predicted))
    predicted_DF = dataMatrix.loc[reads_cov,var_predicted]
    prop_count = pf.compute_proportions(predicted_DF)

    pred_prop_count = pf.create_dictionary(var_predicted, prop_count)
    return pred_prop_count, all_solutions

#A function to run the ADP on individual genes
def run_ADP(editDist, iteration):
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    ADP_dict = defaultdict(int)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE ALLELE DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    ADP_alleles = []
    ADP_prop = []
    for gene in genes:
        print('now processing gene: {}\n'.format(gene))
        result, all_solutions = Process_reads(gene, editDist, iteration)
        print("All the solutions are are: {}\n".format(all_solutions))
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

#run the main function
if __name__ == "__main__":
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
    ap.add_argument("-emd", "--emdDistMode", required=False, default='snp', help="The distance used for EMD calculation. 'snp' compares strains by bases, 'hamming' by genes. Default = 'snp'.")
    ap.add_argument("-art", "--art", required=False, default="art_illumina", help="Path to art_illumina binary file. Default is 'art_illumina'.")
    args = vars(ap.parse_args())

    directories = [d for d in os.listdir(".") if os.path.isdir(d)]
    if args['masterDir'] not in directories:
        os.mkdir(args['masterDir'])

    precision = []
    soft_precision = []
    #weight_list = []
    #numRand_list = []
    emd = []
    soft_recall = []
    recall = []
    total_var_dist = []
    ADP_prec_vals = []
    ADP_rec_vals = []
    ADP_tvd_vals = []
    #known_prec_list = []
    #known_rec_list = []
    if args["experimentType"] == 1:
        if args["modType"] == "Type_1":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, soft_prec, soft_rec, EMD = EvoMod1(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"],args["emdDistMode"], args["art"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                soft_precision.append(soft_prec)
                soft_recall.append(soft_rec)
                emd.append(EMD)
                #weight_list.append(weight)
                #numRand_list.append(numRand)
        elif args["modType"] == "Type_2":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, soft_prec, soft_rec, EMD = EvoMod2(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"], args["emdDistMode"], args["art"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                soft_precision.append(soft_prec)
                soft_recall.append(soft_rec)
                emd.append(EMD)
                #weight_list.append(weight)
                #numRand_list.append(numRand)
        elif args["modType"] == "Type_3":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, soft_prec, soft_rec, EMD = EvoMod3(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"], args["emdDistMode"], args["art"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                soft_precision.append(soft_prec)
                soft_recall.append(soft_rec)
                emd.append(EMD)
                #weight_list.append(weight)
                #numRand_list.append(numRand)

    elif args["experimentType"] == 2:
        if args["modType"] == "Type_1":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, soft_prec, soft_rec, EMD, ADP_prec, ADP_rec, ADP_tvd = EvoMod1(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"], args["emdDistMode"], args["art"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                soft_precision.append(soft_prec)
                soft_recall.append(soft_rec)
                emd.append(EMD)
                #weight_list.append(weight)
                #numRand_list.append(numRand)
                ADP_prec_vals.append(ADP_prec)
                ADP_rec_vals.append(ADP_rec)
                ADP_tvd_vals.append(ADP_tvd)
        elif args["modType"] == "Type_2":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, soft_prec, soft_rec, EMD, ADP_prec, ADP_rec, ADP_tvd, known_prec, known_rec = EvoMod2(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"], args["emdDistMode"], args["art"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                soft_precision.append(soft_prec)
                soft_recall.append(soft_rec)
                emd.append(EMD)
                #weight_list.append(weight)
                #numRand_list.append(numRand)
                ADP_prec_vals.append(ADP_prec)
                ADP_rec_vals.append(ADP_rec)
                ADP_tvd_vals.append(ADP_tvd)
                #known_prec_list.append(known_prec)
                #known_rec_list.append(known_rec)
        elif args["modType"] == "Type_3":
            for i in range(args["numOfiter"]):
                print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Now running Iteration {} based on Experiment {} of Type {}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format((i+1),args["modType"],args["experimentType"])
                pre, rec, tvd, soft_prec, soft_rec, EMD, ADP_prec, ADP_rec, ADP_tvd = EvoMod3(reference,args["editDist"],args["numMut"],i,args["modType"], args["experimentType"], args["emdDistMode"], args["art"])
                precision.append(pre)
                recall.append(rec)
                total_var_dist.append(tvd)
                soft_precision.append(soft_prec)
                soft_recall.append(soft_rec)
                emd.append(EMD)
                #weight_list.append(weight)
                #numRand_list.append(numRand)
                ADP_prec_vals.append(ADP_prec)
                ADP_rec_vals.append(ADP_rec)
                ADP_tvd_vals.append(ADP_tvd)


    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    stats=0
    #print("Statistics {}".format(stats+1))
    #stats += 1
    #avg_kprec = float(sum(known_prec_list))/len(known_prec_list)
    #std_kprec = np.std(np.array(known_prec_list))
    #print 'Precision for known strains:', known_prec_list, "\n"
    #print "Average of precision is: ", avg_kprec, "\n"
    #print "Standard deviation of precision is: ", std_kprec, "\n"

    #print("Statistics {}".format(stats+1))
    #stats += 1
    #avg_krec = float(sum(known_rec_list))/len(known_rec_list)
    #std_krec = np.std(np.array(known_rec_list))
    #print 'Recall for known strains:', known_rec_list, "\n"
    #print "Average of recall is: ", avg_krec, "\n"
    #print "Standard deviation of recall is: ", std_krec, "\n"

    print("Statistics {}".format(stats+1))
    stats += 1
    avg_soft_prec = float(sum(soft_precision))/len(soft_precision)
    std_soft_prec = np.std(np.array(soft_precision))
    print 'Soft Precision is:', soft_precision, "\n"
    print "Average of Soft Precision is: ", avg_soft_prec, "\n"
    print "Standard deviation of Soft Precision is: ", std_soft_prec, "\n"

    print("Statistics {}".format(stats+1))
    stats += 1
    avg_soft_rec = float(sum(soft_recall))/len(soft_recall)
    std_soft_rec = np.std(np.array(soft_recall))
    print 'Soft Recall is:', soft_recall, "\n"
    print "Average of Soft Recall is: ", avg_soft_rec, "\n"
    print "Standard deviation of Soft Recall is: ", std_soft_rec, "\n"

    print("Statistics {}".format(stats+1))
    stats += 1
    avg_emd = float(sum(emd))/len(emd)
    std_emd = np.std(np.array(emd))
    print 'Earth Movers distance is:', emd, "\n"
    print "Average of EMD is: ", avg_emd, "\n"
    print "Standard deviation of EMD is: ", std_emd, "\n"
    
    print("Statistics {}".format(stats+1))
    stats += 1
    avg_prec = sum(precision)/len(precision)
    std_prec = np.std(np.array(precision))
    print 'Precision is:', precision, "\n"
    print "Average of precision is: ", avg_prec, "\n"
    print "Standard deviation of precision is: ", std_prec, "\n"


    print("Statistics {}".format(stats+1))
    stats += 1
    avg_rec = sum(recall)/len(recall)
    std_rec = np.std(np.array(recall))
    print 'Recall is:', recall, "\n"
    print "Average of recall is: ", avg_rec, "\n"
    print "Standard deviation of recall is: ", std_rec, "\n"

    print("Statistics {}".format(stats+1))
    stats += 1
    avg_TVD = sum(total_var_dist)/len(total_var_dist)
    std_TVD = np.std(np.array(total_var_dist))
    print 'Total Variation Distance is:', total_var_dist, "\n"
    print "Average of Total Variation Distance is: ", avg_TVD, "\n"
    print "Standard deviation of Total Variation Distance is: ", std_TVD, "\n"

    if args["experimentType"] == 2:
        os.chdir(args["masterDir"])

        print("Statistics {}".format(stats+1))
        stats += 1
        avg_ADP_prec = sum(ADP_prec_vals)/len(ADP_prec_vals)
        std_ADP_prec = np.std(np.array(ADP_prec_vals))
        print 'Precision for ADP is:', ADP_prec_vals, "\n"
        print "Average of ADP precision is: ", avg_ADP_prec, "\n"
        print "Standard deviation of ADP precision is: ", std_ADP_prec, "\n"

        ADPprecisionDF = pd.DataFrame(ADP_prec_vals)
        ADPprecisionDF = ADPprecisionDF.T
        ADPprecisionDF.to_csv('ADP_{}_editDist_{}_Precision_values.csv'.format(args["modType"], args["editDist"]),sep = '\t')

        print("Statistics {}".format(stats+1))
        stats += 1
        avg_ADP_rec = sum(ADP_rec_vals)/len(ADP_rec_vals)
        std_ADP_rec = np.std(np.array(ADP_rec_vals))
        print 'Recall for ADP is:', ADP_rec_vals, "\n"
        print "Average of ADP recall is: ", avg_ADP_rec, "\n"
        print "Standard deviation of ADP recall is: ", std_ADP_prec, "\n"

        ADPrecallDF = pd.DataFrame(ADP_rec_vals)
        ADPrecallDF = ADPrecallDF.T
        ADPrecallDF.to_csv('ADP_{}_editDist_{}_Recall_values.csv'.format(args["modType"],args["editDist"]), sep = '\t')

        print("Statistics {}".format(stats+1))
        stats += 1
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


    #EMD
    plt.figure()
    plt.hist(emd, bins=np.linspace(-1,1))
    plt.title("Plot of EMD for edit distance {}".format(args["modType"], args["editDist"]))
    plt.xlabel("EMD")
    plt.ylabel("Frequency")
    plt.savefig('{}_editDist_{}_EMD_plot'.format(args["modType"],args["editDist"]))
    #Save the plot for boxplot plotting
    emdDF = pd.DataFrame(emd)
    emdDF = emdDF.T
    emdDF.to_csv('{}_editDist_{}_EMD_values.csv'.format(args["modType"],args["editDist"]), sep = '\t')

    #Soft recall
    plt.figure()
    plt.hist(soft_recall, bins=np.linspace(0,2))
    plt.title("Plot of simulation Soft Recall for {} for edit Distance {}".format(args["modType"], args["editDist"]))
    plt.xlabel("Soft Recall")
    plt.ylabel("Frequency")
    plt.savefig('{}_editDist_{}_SoftRecall_plot'.format(args["modType"],args["editDist"]))
    #Save the plot for boxplot plotting
    softrecallDF = pd.DataFrame(soft_recall)
    softrecallDF = softrecallDF.T
    softrecallDF.to_csv('{}_editDist_{}_SoftRecall_values.csv'.format(args["modType"],args["editDist"]), sep = '\t')

    #Soft precision
    plt.figure()
    plt.hist(soft_precision, bins=np.linspace(0,2))
    plt.title("Plot of simulation Soft Precision for {} for edit Distance {}".format(args["modType"], args["editDist"]))
    plt.xlabel("Soft Precision")
    plt.ylabel("Frequency")
    plt.savefig('{}_editDist_{}_SoftPrecision_plot'.format(args["modType"], args["editDist"]))
    #save the plot for boxplot plotting
    soft_precisionDF = pd.DataFrame(soft_precision)
    soft_precisionDF = soft_precisionDF.T
    soft_precisionDF.to_csv('{}_editDist_{}_SoftPrecision_values.csv'.format(args["modType"], args["editDist"]),sep = '\t')

    #recall
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
