#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 16:14:39 2017

@author: elijah
This scripts similuate a set of strains randomly chosen with random proportions as well.
We are also looking to introduce new strains as errors
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
plt.style.use('ggplot')


seed = 1994
numpy.random.seed(seed)

def Mutate_strain(reference,hamming_dist,loci,Type):
    '''
        A function that takes a given strain and mutates it according to a given
        hamming distance

        reference: The reference database table of all the currrently known strains

        hamming_dist: The hamming distance/number of mutations to subject a given strain to

        loci: A list containing all the loci
    '''
    final_strains = []
    #pick a random strain to mutate
    if Type == 1:
        k = random.randint(1,reference.shape[0])
        strain = reference.iloc[k,:].tolist()
        final_strains.append(strain)
        #pick random positions on the strain to mutate
	random_index = np.random.choice(7,hamming_dist,replace=False)
        #random_index = random.sample(xrange(0,7),hamming_dist)
        temp = defaultdict(list)
        #list to hold the new strain
        new_strain = [None]*8
        #generate a random strain by selecting a random allele for a selected gene from the database and replacing it
        for index in random_index:
            locus = reference[loci[index]].tolist()
            locus = list(set(locus))
            locus.remove(strain[index])
            new_allele = random.choice(locus)
            temp[index] = new_allele
        #update the strains
        for i in range(len(loci)):
            if i in random_index:
                new_strain[i] = temp[i]
            else:
                new_strain[i] = strain[i]
        final_strains.append(new_strain)
    else:
        strains = [] #list to hold the strains
	strain_index = np.random.choice(729,3,replace=False)+1
        #strain_index = random.sample(xrange(1,reference.shape[0]),3) #randomly select which three strains to use
        for index in strain_index: #select the strains
            strains.append(reference.iloc[index,:].tolist())
            final_strains.append(reference.iloc[index,:].tolist())
        random_strains = random.sample(strains,2) #randomly select two strains to mutate
        for rand_strain in random_strains:
            temp = defaultdict(list)
	    random_index = np.random.choice(7,hamming_dist,replace=False)
            #random_index = np.random.random_integers(0,7,hamming_dist)
            new_strain = [None]*8
            for index in random_index:
                locus = reference[loci[index]].tolist()
                locus = list(set(locus))
                locus.remove(rand_strain[index])
                new_allele = random.choice(locus)
                temp[index] = new_allele
            for i in range(len(loci)):
                if i in random_index:
                    new_strain[i] = temp[i]
                else:
                    new_strain[i] = rand_strain[i]
            final_strains.append(new_strain)


    #return both the random strain and its mutated counterpart
    return final_strains

def check_strain(strain,strains):
    if strain not in strains:
        return True
def create_strain(strain1,strain2):
    new_strain = []
    for i in range(1,9):
    	val = np.random.choice(np.arange(1, 3), p=[0.5,0.5])
        if val == 1:
            new_strain.append(strain1[i-1])
        else:
            new_strain.append(strain2[i-1])
    return new_strain

def Recomb_strains(reference,Type):
    '''
        This function takes the reference database, picks two strains and randomly
        recombine them both to create and new and third strain not in the database.

        Reference: The reference database table of all the currrently known strains
    '''
    final_strains = []
    if Type == 1:
        #randomly pick two strains from the reference database
	randomStrainIndex= np.random.choice(729,2,replace=False)+1
        #randomStrainIndex = random.sample(xrange(1,731),2)
        strains = []
        new_strain = []
        #get the new strains
        for num in range(len(randomStrainIndex)):
            strains.append(reference.iloc[randomStrainIndex[num],:].tolist())
            final_strains.append(reference.iloc[randomStrainIndex[num],:].tolist())
        #randomly choose which allele to include in the new strain with equal probability
        for i in range(1,9):
            val = np.random.choice(np.arange(1, 3), p=[0.5,0.5])
            if val == 1:
                new_strain.append(strains[0][i-1])
            else:
                new_strain.append(strains[1][i-1])
        final_strains.append(new_strain)

    else:
        #pick three random strains to mutate
        strains = [] #list to hold the strains
	strain_index = np.random.choice(729,3,replace=False)+1
        #strain_index = random.sample(xrange(1,reference.shape[0]),3) #randomly select which three strains to use
        for index in strain_index: #select the strains
            strains.append(reference.iloc[index,:].tolist())
            final_strains.append(reference.iloc[index,:].tolist())
        #generate all pairs for recombination
        strain_pairs = [pair for pair in itertools.combinations(strains,r=2)]
        for pair in strain_pairs:
            new_strain = []
            temp = list(pair)
            '''for i in range(1,9):
                val = np.random.choice(np.arange(1,3), p=[0.5,0.5])
                if val == 1:
                    new_strain.append(temp[0][i-1])
                else:
                    new_strain.append(temp[1][i-1])
     	    final_strains.append(new_strain)'''
            flag = True
	    new_strain = create_strain(temp[0],temp[1])
	    while flag:
                if check_strain(new_strain,final_strains):
		    final_strains.append(new_strain)
		    flag = False
 		else:
 		    new_strain = create_strain(temp[0],temp[1])
                    flag = True
    return final_strains


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
    print type(strain_df)
    for i in range(strain_df.shape[0]):
        row = strain_df.iloc[i].tolist()
        temp_key = tuple(row[0:8])
        predicted_dict[temp_key] = float(row[10])
    print ('The True strains are {}'.format(strain_dict.keys()))
    print ('The True proportions are {}'.format(strain_dict.values()))
    print ('The Predicted strains are {}'.format(predicted_dict.keys()))
    print ('The Predicted proportions are {}'.format(predicted_dict.values()))
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

def Simulate_strains(numIter,strainRef,sampleDir,outputDir,samplesDir,simFilesDir,sim_type,hamming_dist,Type):
    '''
        numIter: The number of simulation iterations

        strainRef: The reference database table of all the currrently known strains

        sampleDir: The directory where the temporary files for the current sample will be stored

        outputDir: The Directory where the outputs will be stored

        samplesDIr: The directory containing all the samples

        sim_type: The type of simulation being run

        hamming_dist: The number of mutations to introduce into the new strain

        simFilesDir: The directory to store the intermediate files generated
    '''

    print ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~You have chosen to run {0} case simulation for type {1}. Now running simulation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'.format(Type,sim_type))
    if sim_type == 'Mutation':
        if Type == 1:
            numStrains = 2
        else:
            numStrains = 5
    if sim_type == 'Recombination':
        if Type == 1:
            numStrains = 3
        else:
            numStrains = 6
    #loci used to generate strains
    loci = ['clpA','clpX','nifS','pepX','pyrG','recG','rplB','uvrA']
    precision = []
    recall = []
    total_var_dist = []
    times = []
    #read in the current know strain database
    reference = pd.read_csv(strainRef,sep="\t",usecols=range(1,len(loci)+1))
    #parse the database
    for name in loci:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
    for iteration in range(1,numIter+1):
        #start the timer
        start = time.time()
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running Simulation {} ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~".format(iteration))
        #Go to the current sample directory
        os.chdir(sampleDir)

        #define dictionaries to hold each of the variants being used in the simulated strains
        clpA = defaultdict(list)
        clpX = defaultdict(list)
        nifS = defaultdict(list)
        pepX = defaultdict(list)
        pyrG = defaultdict(list)
        recG = defaultdict(list)
        rplB = defaultdict(list)
        uvrA = defaultdict(list)
	aux_strain_index =  np.random.choice(729,1,replace=False)+1
	aux_strain = reference.iloc[aux_strain_index,:].values.tolist()
	aux_strain = tuple(aux_strain[0])
        #list to hold the names of the varaint dictionaries
        dict_list = [clpA,clpX,nifS,pepX,pyrG,recG,rplB,uvrA]
        #dictionary to hold the strains randomly simulated
        strain_dict = defaultdict(list)
        #generate random proportions for the strains
        nums = [random.uniform(0,1) for x in range(0,numStrains)]
        Sum = reduce(lambda x,y: x+y, nums)
        proportions = [x/Sum for x in nums]

	#proportions = np.random.dirichlet(np.ones(numStrains),1)
	#print sum(proportions[0])
        print sum(proportions)
        #proportions = [random.random() for j in range(numStrains)]
        #prop_sum = sum(proportions)
        #proportions = [i/prop_sum for i in proportions]
        #define a dictionary of dictionaries to hold all the strains and their respective proportions.
        strain_var_proportions = defaultdict(dict)
        #mutate the strain to generate a new strain
        if sim_type == 'Mutation':
            strain_list = Mutate_strain(reference,hamming_dist,loci,Type)
        #recombine two strains to generate a third(new) one
        if sim_type == 'Recombination':
            strain_list = Recomb_strains(reference,Type)
        #update the proportions of each gene allele in the strains
        for i in range(len(strain_list)):
            strain_dict[tuple(strain_list[i])] = proportions[i]#*100
            strain_var_proportions[tuple(strain_list[i])] = Compute_Variant_proportions(strain_list[i],proportions[i])
        #update the gene dictionaries with the appropraite proportions
        for strain_key in strain_var_proportions.keys():
            #define a dictionary to hold the current strain and its proportions
            gene_dict = strain_var_proportions[strain_key]
            for gene_key in gene_dict:
                #for clpA
                if 'clpA' in gene_key:
                    if gene_key not in clpA:
                        clpA[gene_key] = gene_dict[gene_key]
                    else:
                        clpA[gene_key] += gene_dict[gene_key]
                #for clpX
                if 'clpX' in gene_key:
                    if gene_key not in clpX:
                        clpX[gene_key] = gene_dict[gene_key]
                    else:
                        clpX[gene_key] += gene_dict[gene_key]
                #for nifS
                if 'nifS' in gene_key:
                    if gene_key not in nifS:
                        nifS[gene_key] = gene_dict[gene_key]
                    else:
                        nifS[gene_key] += gene_dict[gene_key]
                #for pepX
                if 'pepX' in gene_key:
                    if gene_key not in pepX:
                        pepX[gene_key] = gene_dict[gene_key]
                    else:
                        pepX[gene_key] += gene_dict[gene_key]
                #for pyrG
                if 'pyrG' in gene_key:
                    if gene_key not in pyrG:
                        pyrG[gene_key] = gene_dict[gene_key]
                    else:
                        pyrG[gene_key] += gene_dict[gene_key]
                #for recG
                if 'recG' in gene_key:
                    if gene_key not in recG:
                        recG[gene_key] = gene_dict[gene_key]
                    else:
                        recG[gene_key] += gene_dict[gene_key]
                #for rplB
                if 'rplB' in gene_key:
                    if gene_key not in rplB:
                        rplB[gene_key] = gene_dict[gene_key]
                    else:
                        rplB[gene_key] += gene_dict[gene_key]
                #for uvrA
                if 'uvrA' in gene_key:
                    if gene_key not in uvrA:
                        uvrA[gene_key] = gene_dict[gene_key]
                    else:
                        uvrA[gene_key] += gene_dict[gene_key]

	print sum(clpA.values())
	print sum(clpX.values())
	print sum(nifS.values())
	print sum(pepX.values())
	print sum(pyrG.values())
	print sum(recG.values())
	print sum(rplB.values())
	print sum(uvrA.values())
        #if sum(clpA.values()) != 1:
            #continue
        filesHere = [d for d in os.listdir(".") if os.path.isfile(d)]
        for csv_file in filesHere:
            os.rename(csv_file,simFilesDir+'/'+csv_file)
        #write the dictionaries to a csv file
        for i in range(len(dict_list)):
            locus = dict_list[i]
            Dict_to_csv(locus,loci[i])

        #Compute the simulation statistics
        os.chdir(samplesDir)
        #os.remove(samplesDir+'/borreliaLP.lp')
        strain_df = pf.strainSolver(samplesDir,strainRef,outputDir,loci,'all','all',10800,5)
        pre,rec,tvd = Compute_Prec_and_rec(strain_dict,strain_df)
        #compute the time taken
        time_taken = time.time() - start
        times.append(time_taken)
        precision.append(pre)
        recall.append(rec)
        total_var_dist.append(tvd)
        os.chdir(sampleDir)
        filesHere = [d for d in os.listdir(".") if os.path.isfile(d)]
        for csv_file in filesHere:
            os.rename(csv_file,simFilesDir+'/'+csv_file)

    return precision,recall,total_var_dist,times



def main():

    #get and parse the arguments
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--numOfiter", required=False, default=40, type=int, help="The number of simulation iterations default = 40")
    ap.add_argument("-d", "--masterDir", required=False, default='Results',type=str, help="The name of the master directory to run the simulation and store the results. Default = Created folder called Results")
    ap.add_argument("-r", "--strainRef", required=True, help="Name of the text file containing the reference strains that is stored in this directory")
    ap.add_argument("-t", "--Sim_Type", required=False, default='Mutation', type=str, help="The type of simulation you would like to run (Mutation or Recombination). Default = Mutation")
    ap.add_argument("-hd", "--HammingDist", required=False, default=2, type=int, help="The total number of mutations you would like to make to a strain. Default = 2")
    ap.add_argument("-st", "--Mut_Rec_Type", required=False, default=1, type=int, help="The type of simulation to run. Default = 1 which is the simple case, and 2 = the complex case")
    args = vars(ap.parse_args())

    if args["Sim_Type"] == "Mutation":
        outputDir = "{0}_{1}_Simulation_Output".format(args['Sim_Type'],args["HammingDist"])
        sampleDir = "{0}_{1}_Simulation_Samples".format(args['Sim_Type'],args["HammingDist"])
        sample = "{0}_{1}_Simulation_001".format(args['Sim_Type'],args["HammingDist"])
        tempDir = "{0}_{1}_Simulation_files".format(args["Sim_Type"],args["HammingDist"])
    elif args["Sim_Type"] == "Recombination":
        outputDir = "{}_Simulation_Output".format(args['Sim_Type'])
        sampleDir = "{}_Simulation_Samples".format(args['Sim_Type'])
        sample = "{}_Simulation_001".format(args['Sim_Type'])
        tempDir = "{}_Simulation_files".format(args["Sim_Type"])
    directories = [d for d in os.listdir(".") if os.path.isdir(d)]
    strainRef = os.path.join(str(os.getcwd()),"strain_ref.txt")
    if args['masterDir'] not in directories:
        os.mkdir(args['masterDir'])
    masterDir = os.path.join(str(os.getcwd()), args['masterDir'])
    os.chdir(args["masterDir"]) #change into the master directory
    directoriesHere = [d for d in os.listdir(".") if os.path.isdir(d)]
    if sampleDir not in directoriesHere and sample not in directoriesHere and outputDir not in directoriesHere and tempDir not in directoriesHere:
        os.mkdir(sampleDir) #make a directory containing all the samples
        os.mkdir(outputDir) #Make a directory for the outputs
        os.mkdir(tempDir) #directory to hold all the simulation files
        os.chdir(sampleDir)
        os.mkdir(sample) #make a directory for the current simulation sample
    else:
        os.chdir(tempDir)
        filesHere = [d for d in os.listdir(".") if os.path.isfile(d)]
        for old_file in filesHere:
            os.remove(old_file)
    
    #run the simulation
    precision,recall,total_var_dist,times = Simulate_strains(args["numOfiter"],strainRef,os.path.join(masterDir, sampleDir, sample),os.path.join(masterDir, outputDir),os.path.join(masterDir, sampleDir),os.path.join(masterDir, tempDir),args['Sim_Type'],args["HammingDist"],args["Mut_Rec_Type"])
    my_list = [precision,recall,total_var_dist]
    my_df = pd.DataFrame(my_list)
    my_df.to_csv('Values_for_plot.csv')
    #plot the results
    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Creating plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    
    os.chdir(os.path.join(masterDir, outputDir))
    plt.figure()
    plt.hist(recall, bins=np.linspace(0,2))
    plt.title("Plot of simulation recall for {}".format(args['Sim_Type']))
    plt.xlabel("Recall")
    plt.ylabel("Frequency")
    plt.savefig('Recall_plot')
    #Save the plot for boxplot plotting
    recallDF = pd.DataFrame(recall)
    recallDF = recallDF.T
    recallDF.to_csv('Recall_values.csv', sep = '\t')
    #for precision
    plt.figure()
    plt.hist(precision, bins=np.linspace(0,2))
    plt.title("Plot of simulation precision for {}".format(args['Sim_Type']))
    plt.xlabel("Precision")
    plt.ylabel("Frequency")
    plt.savefig('Precision_plot')
    #save the plot for boxplot plotting
    precisionDF = pd.DataFrame(precision)
    precisionDF = precisionDF.T
    precisionDF.to_csv('Precision_values.csv',sep = '\t')
    #for total variation distance
    plt.figure()
    plt.hist(total_var_dist,bins=np.linspace(-1,1))
    plt.title("Plot of simulation Total variation distance for {}".format(args['Sim_Type']))
    plt.xlabel("Total Variation Distance")
    plt.ylabel("Frequency")
    plt.savefig('Total_variation_Distance_plot')
    #save the plot for boxplot plotting
    total_var_distDF = pd.DataFrame(total_var_dist)
    total_var_distDF = total_var_distDF.T
    total_var_distDF.to_csv('Total_Variation_Distance.csv', sep = '\t')
    #for the time taken
    plt.figure()
    plt.hist(times)
    plt.title("Plot of simulation Total time taken for {}".format(args['Sim_Type']))
    plt.xlabel("Total time taken")
    plt.ylabel("Frequency")
    plt.savefig('Totat_time_taken_plot')

    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

    print ("The recall values for all simulations are: {}".format(recall))
    print ("The precision values for all simulations are: {}".format(precision))
    print ("The total variation distances for all simulations are: {}".format(total_var_dist))
    print ("The time taken for the simulations are: {}".format(times))

#run the main function
if __name__ == "__main__":
    main()
