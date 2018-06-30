#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  6 15:51:11 2018

@author: elijah
"""

from __future__ import division
from collections import defaultdict
import pipeline_functions as pf
import variantILP as varSolver
import pandas as pd
import os
import random
import sh
import csv
import sys
import numpy as np

#Global variable
_EDITDIST = 5    #For soft precision and recall calculation
MAX_WEIGHT = 100

#get get the current working directory
curr_dir = os.getcwd()
    
#get the path to the loci dist db
dist_db = os.path.join(curr_dir, 'Dist_DB')

#set up genes edit distance matrices
#clpA
clpA_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_clpA.csv'), sep=",")
#clpX
clpX_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_clpX.csv'), sep=",")
#nifS
nifS_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_nifS.csv'), sep=",")
#pepX
pepX_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_pepX.csv'), sep=",")
#pyrG
pyrG_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_pyrG.csv'), sep=",")
#recG
recG_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_recG.csv'), sep=",")
#rplB
rplB_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_rplB.csv'), sep=",")
#uvrA
uvrA_df = pd.read_csv(os.path.join(dist_db, 'editDistanceMatrix_uvrA.csv'), sep=",")

genes_df = [clpA_df,clpX_df,nifS_df,pepX_df,pyrG_df,recG_df,rplB_df,uvrA_df]

loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']

genes_dict = {locus:[genes_df[i]] for i,locus in enumerate(loci)}


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


#function that mutates a strain given an edit dist and number of mutations
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
    
#check if a strain is already in the database    
def check_strain(strain,strains):
    if strain not in strains:
        return True
    else:
        return False
    

def Recombine_strains(strain1, strain2):
    '''
    A function that takes in two strains and return a recombination betweent the two
    
    Input: Strain1 and Strain2
    
    Output: Strain3, a recombination between the two strains
    '''
    result = []
    for i in range(1,9):
        pos = np.random.choice(np.arange(1, 3), p=[0.5,0.5])
        if pos == 1:
            result.append(strain1[i-1])
        else:
            result.append(strain2[i-1])
    return result        

        
def Generate_reads(strains, seed, editDist):
    variants_path = os.path.join(curr_dir, 'Variant_files')
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    nums = [random.uniform(0,1) for x in range(0,len(strains))]
    Sum = reduce(lambda x,y: x+y, nums)
    proportions = [x/Sum for x in nums]
    for i in range(len(strains)):
        sequence_file = open('editDist_{}_Strain_{}_sequences.fas'.format(editDist,i), 'w')
        strain = strains[i]
        output_file_name = 'editDist_{}_strain_{}_reads_'.format(editDist,i)
        for j in range(len(strain)):
            variant_sequence = sh.grep(str(strain[j]),"{}/{}_linear.txt".format(variants_path, genes[j]),"-w","-A1") #use bash to extract the variant sequence
            variant_sequence = variant_sequence.rstrip() #remove the "\n" character that is returned by bash
            variant_sequence = str(variant_sequence)
            #total_variants_sequence += variant_sequence
            sequence_file.write('{0}{1}'.format(variant_sequence, '\n'))
        #run art to generate reads with 100X coverage
        #store
        sequence_file.close()
        #Set the ART command, I have included a random seed for reproducibility, and a coverage parameter
        ART_command = "art_illumina -qL 33 -qs 10 -qs2 15 -k 3 -rs {0} -q -ss HS25 -sam -i {1} -p -l 76 -f {2} -m 149.452058562 -s 41.3519698623 -o {3}".format(seed, 'editDist_{}_Strain_{}_sequences.fas >/dev/null 2>&1'.format(editDist,i),proportions[i]*100, output_file_name)
        os.system(ART_command)
    appendFirst_cmd = "cat *_1.fq > editDist_{}_all_Strains_1.fa".format(editDist) #append all the first of the pairs together
    appendSecond_cmd ="cat *_2.fq > editDist_{}_all_Strains_2.fa".format(editDist) #append all the second of the pairs together
    os.system(appendFirst_cmd)
    os.system(appendSecond_cmd)
    os.system('rm *aln*')
    os.system('rm *seq*')
    os.system('rm *sam*')
    #os.system('rm *fa*')
    os.system('rm *fq*')
    return proportions          
        

    
#A function to run the ADP on individual genes
def run_ADP(editDist):
    genes = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
    ADP_dict = defaultdict(int)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE ALLELE DETECTION ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    ADP_alleles = []
    ADP_prop = []
    for gene in genes:
        print('now processing gene: {}\n'.format(gene))
        result, all_solutions = Process_reads(gene, editDist)
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
        #print('The sum of the alleles proportionsss are: {}'.format(sum(result.values())))
        #os.chdir('../')
        #ADP_dict[gene] = result
    return ADP_dict, ADP_alleles, ADP_prop
    
    
    
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
    
            
            
def Process_reads(gene, editDist):
    loci_db = os.path.join(curr_dir, 'Bowtie_Indices')
    mapping_cmd = "bowtie -a --best --strata -v 3 -p 4 {0}/{1}_bowtie -1 ./editDist_{2}_all_Strains_1.fa -2 ./editDist_{2}_all_Strains_2.fa --sam all_Strains.sam".format(loci_db, gene, editDist)
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
    pred_object_val,var_predicted,reads_cov,all_solutions, all_objective = varSolver.solver(dataMatrix,Qmatrix, None)
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
            tvd = tvd + abs(true[key] - predicted[key])
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
    
def Compute_avg(precision,recall, tvd):
    
    avg_prec = sum(map(float, [item for sublist in precision.values() for item in sublist]))/len(precision.values())
    avg_rec = sum(map(float, [item for sublist in recall.values() for item in sublist]))/len(precision.values())
    avg_tvd = sum(map(float, [item for sublist in tvd.values() for item in sublist]))/len(precision.values())
    
    return avg_prec, avg_rec, avg_tvd
    

def Write_Proportions(strains, proportions):
    #create a strain and its proportions dict
    strain_prop_dict = defaultdict(list)
    for i in range(len(strains)):
        #print tuple(strains[i])
        strain_prop_dict[tuple(strains[i])] = proportions[i]
    return strain_prop_dict

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
        temp_key = tuple(row[2:10])
        #print(row[9])
        predicted_dict[temp_key] = float(row[10])
    #print ('The True strains and their proportions are {}\n'.format(strain_dict))
    #print ('The Predicted strains and their proportions are\n {}'.format(predicted_dict))
    #Compute the total variation distance
    for key in predicted_dict.keys():
        if key in strain_dict.keys():
            total_var_dist += abs(round(strain_dict[key],3)-predicted_dict[key])
            count += 1
        else:
            total_var_dist += predicted_dict[key]
    for key in strain_dict.keys():
        if key not in predicted_dict.keys():
            #print(strain_dict[key])
            total_var_dist += strain_dict[key]
    #compute the precision
    precision = count/strain_df.shape[0]
    #compute the recall
    recall = count/len(strain_dict)

    return precision, recall, total_var_dist/2.0

def Compute_SEST_Soft_Prec_and_Rec(true, predicted, editDist):
    PredictedStrains = [ tuple(sorted(i)) for i in predicted.keys() ]
    TrueStrains = [ tuple(sorted(i)) for i in true.keys() ]
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
    
    
def EvoMod1(reference, editDist, num_mut, iteration, all_strains, seed, art_cmd="art_illumina", pathToDistMat=None):
    strains = []
    #select random strain from the first 100
    S1_index =  random.randint(1,99)
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
            S2 = Mutate_strain(S1, editDist, num_mut)
            flag = True
    #get the true alleles
    true_alleles = set([item for sublist in strains for item in sublist])
    true_all = defaultdict(int)

    #Generate the reads
    proportions = Generate_reads(strains, seed, editDist)
    for strain in strains:
        strain = list(strain)
        for gene in strain:
            if gene not in true_all:
                true_all[gene] = proportions[strains.index(strain)]
            else:
                true_all[gene] = true_all[gene] + proportions[strains.index(strain)]
    return proportions, true_all, true_alleles, strains
    

def EvoMod2(reference, editDist, num_mut, iteration, all_strains, seed, art_cmd="art_illumina", pathToDistMat=None, dropStrain=None):
    strains = []
     #get the indices of the strains to use
    indices = random.sample(range(1, 98), 2)
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
                S2 = Mutate_strain(S1, editDist, num_mut)
                flag = True
    true_all = defaultdict(int)
    
    #EvoMod2.1: Drop an existing strain, resulting in 3 strains
    if dropStrain == "existing":
        existingInd = [i for i in range(len(strains)) if i % 2 == 0 ]
        removeInd = random.sample(existingInd, 1)[0]
        strains = [st for idx, st in enumerate(strains) if idx != removeInd]
    elif dropStrain == "new":
        newInd = [i for i in range(len(strains)) if i % 2 == 1 ]
        removeInd = random.sample(newInd, 1)[0]
        strains = [st for idx, st in enumerate(strains) if idx != removeInd]
    #Generate the reads
    proportions = Generate_reads(strains, seed, editDist)
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
    true_alleles = set([item for sublist in strains for item in sublist])
    return proportions, true_all, true_alleles, strains
    
def EvoMod3(reference, editDist, num_mut, iteration, all_strains, seed, art_cmd="art_illumina", pathToDistMat=None):
    '''
    EvoMod3(Result in 5 strains, 3 new 2 existing)
    1) Do EvoMod2
    2) Choose any 2 strains to recombine, to produce a 3rd new strain
    3) If we can afford to do, 2 recombination events. We could try that
    4) ...
    '''
    
    #do EveMod2
    strains = []
    #get the indices of the strains to use
    indices = random.sample(range(1, 99), 2)
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
                S2 = Mutate_strain(S1, editDist, num_mut)
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
    proportions = Generate_reads(strains, seed, editDist)
    for strain in strains:
        strain = list(strain)
        for gene in strain:
            if gene not in true_all:
                true_all[gene] = proportions[strains.index(strain)]
            else:
                true_all[gene] = true_all[gene] + proportions[strains.index(strain)]  
    #get the set of true alleles  
    true_alleles = set([item for sublist in strains for item in sublist])
    
    return proportions, true_all, true_alleles, strains
    

def Find_Alleles(file):
    '''
    a function that takes a genome fasta file, extracts and returns the names of the alleles found in that genome  
    
    Input: Path to genome fasta file
    
    Output: List of names of alleles found in said genome
    
    '''
    
    genome = open(file, 'r')
    alleles = []
    for line in genome.readlines():
        if '>' in line:
            line = line.replace('>', '')
            alleles.append(line.rstrip())
    return alleles    
    
    
def Compute_SEST_Prec_and_rec(true, predicted):
    '''
    A function that computes the precision ad recall for the StrainEST alogirthm
    
    Inputs:
        true-A dictionary containing the true strains and their proportions
        predicted-A dictionary containing the predicted strains by StrainEST and their proportions
    Outputs:
        The precision, recall values and total variation distance values
    '''
    count = 0
    total_var_dict = 0
    for strain in true.keys():
        if strain in predicted.keys():
            #update the count if we found a matching pair
            count += 1
            #update the total variation distance
            total_var_dict += abs(true[strain] - predicted[strain])
        else:
            total_var_dict += abs(true[strain])
    for strain in predicted.keys():
        if strain not in true.keys():
            total_var_dict += abs(predicted[strain])
            
    #compute the precision
    prec = count/(len(predicted))
    #compute the recall
    rec = count/(len(true))
    
    return prec, rec, total_var_dict/2.0
    
def Compute_SEST_Allele_prec_and_rec(true, predicted):
    '''
    A function that computes the precision ad recall for the StrainEST alogirthm
    
    Inputs:
        true-A list containing the true alleles
        predicted-A list containing the predicted alleles by StrainEST
    Outputs:
        The precision and recall values
    '''
    count = 0
    
    for allele in true:
        if allele in predicted:
            #update the count if we found a matching pair
            count = count +1
    #compute the precision
    prec = count/(len(predicted))
    #compute the recall
    rec = count/(len(true))
    
    return prec, rec
    

def Run_strainEST(editDist):
    '''
    A function that runs the StrainEST algorithm on simulated data. 
    
    Input: the edit distance being used to simulate the data
    
    Output: The set of alleles which constitute the strains the algorithm found
    '''
    
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~NOW RUNNING THE STRAINEST ALGORITHM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    genomes_dir = os.path.join(curr_dir, 'Variant_files')
    output_dir = os.path.join(curr_dir, 'Analysis')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    database = os.path.join(curr_dir, 'Loci_DB_bowtie2')
    print('Now running the first command')
    cmd1 = 'bowtie2 --very-fast --no-unal -x {0}/MA -1 editDist_{1}_all_Strains_1.fa -2 editDist_{1}_all_Strains_2.fa -S editDist_{1}_all_Strains.sam'.format(database, editDist)
    os.system(cmd1)
    #second command
    print ('Now running the second command')
    cmd2 = 'samtools view -b editDist_{0}_all_Strains.sam > editDist_{0}_all_Strains.bam'.format(editDist)
    os.system(cmd2)
    #third command
    print ('now running the third command')
    cmd3 = 'samtools sort editDist_{0}_all_Strains.bam -o editDist_{0}_all_Strains.Sorted.bam'.format(editDist)
    os.system(cmd3)
    #fourth command
    'Now running the fourth command'
    cmd4 = 'samtools index editDist_{0}_all_Strains.Sorted.bam'.format(editDist)
    os.system(cmd4)
    #fifth command
    'Now running the fifth command'
    cmd5 = 'strainest est {0}/snp_clust.dgrp editDist_{1}_all_Strains.Sorted.bam Analysis'.format(database, editDist)
    os.system(cmd5)
    #parse the results
    results = pd.read_csv(os.path.join(output_dir, 'abund.txt'), header=0, sep='\t')
    #Get the representative genome
    genomes = results.loc[results['editDist_{}_all_Strains.Sorted.bam'.format(editDist)] > 0]
    result_list = []
    #compute the strains in the genome
    #create a dictionary to hold the resulting strains and their proportions
    result_dict = defaultdict(int)
    #create a dict to hold the alleles and their proportions
    allele_dict = defaultdict(int) 

    for genome in genomes.iterrows():
        #extract the genome file
        genome_files = genome[1][0]
        #compute the alleles in the genome
        result_alleles = Find_Alleles(os.path.join(genomes_dir,'{}'.format(str(genome_files))))
        #append all alleles and their proportions to a dictionary
        for allele in result_alleles:
            if allele not in allele_dict.keys():
                allele_dict[allele] = genome[1][1]
            else:
                allele_dict[allele] += genome[1][1]

        #append it to the result list
        result_list.append(result_alleles)
        #update the dictionary with the alleles and its proportions
        result_dict[tuple(result_alleles)] = genome[1][1]
    #compute the alleles
    alleles = set([item for sublist in result_list for item in sublist])
    #return the strains and alleles
    return result_list, alleles, result_dict, allele_dict
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
