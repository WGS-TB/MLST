#!/usr/bin/python
from __future__ import division
import sh
import math
import numpy as np
import random
import os
import itertools
import sys
import norm_variantILP as varSolver
import linecache
import pipeline_functions as pf

#Testing purposes and global variables
TEST_EMPTY_LIST = True

''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Function Definitions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''

'''
Input: pred_prop and true_prop which are dictionaries
Return: Total variation distance
'''
def totalVariationDist(pred_prop, true_prop):
    common_variants = set.intersection(set(pred_prop.keys()), set(true_prop.keys()))
    diff_vars_pred = set.difference(set(pred_prop.keys()), common_variants)
    diff_vars_true = set.difference(set(true_prop.keys()), common_variants)
    common_variants = list(common_variants)
    diff_vars_pred = list(diff_vars_pred)
    diff_vars_true = list(diff_vars_true)
    
    totalVarDist=0
    for var in common_variants:
        totalVarDist += abs(pred_prop[var] - true_prop[var])
        
    for var in diff_vars_pred:
        totalVarDist += pred_prop[var]
        
    for var in diff_vars_true:
        totalVarDist += true_prop[var]
        
    return totalVarDist/2

#Return x with first letter being capital
def upperfirst(x):
    return x[0].upper() + x[1:]

#predicted and true are lists
def precision(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(predicted))
       
#predicted and true are lists 
def recall(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(true))

#Return boolean whether predicted matches true
def predictedCorrectly(predicted, true):
        return set(predicted) == set(true)

def writeReadTable(capGene, iteration, option):
    readOutFile = open("{0}_{1}_{2}NoHeader.sam".format(capGene, iteration, option))
    writefile = open("{0}_{1}_{2}_reads.txt".format(capGene, iteration, option), "w")
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

'''
Input: Dataframe with rows=reads, columns=variants
Output: The objective value of these variants, reads which do not map to true variants if any
'''
def compute_true_objective_val(dataframe):
        objective_val = 0
        #To record reads which do not map to true variants
        bad_read = list()
        #Even no limit on mm, there are reads which are covered by predicted variants but not
        #true variants. Identify largest mm , if this case happens, increment by max_mm+1
        max_mm = max(dataframe.max())
        for row in dataframe.itertuples():
            mmInfo_list = [i for i in list(row)[1:] if i >= 0]
            
            #mmInfo_list will be empty if a read does not map back to the true variants
            if TEST_EMPTY_LIST == True:
                if len(mmInfo_list) > 0:
                    objective_val += min(mmInfo_list)   #Increment by the minimum number of mismatches
                else:
                    objective_val+=max_mm + 1
                    bad_read.append(list(row)[0])
#                    print(list(row))
            else:
                if len(mmInfo_list) == 0:
                    bad_read.append(list(row)[0])
#                    print(list(row))
                    
                objective_val += min(mmInfo_list)
        objective_val += len(dataframe.columns)     #Increment by the number of variants used
        return objective_val, bad_read

#Compute the difference between predicted and true objective values
def compute_obj_diff(predicted, true):
        diff = predicted - true
        return diff

#def outputDataForML(dataMatrix, csvfile, varProp_dict):
#    mm = range(7)
#    matrixForML = list()
#    variants = varProp_dict.keys()
#    print(variants)
#    proportions = varProp_dict.values()
#    print(proportions)
#    
#    track=0
#    for v in variants:
#        temp_array = list()
#        temp = dataMatrix.loc[:, v].value_counts()
#        
#        for i in mm:
#            if i in temp.index:
#                temp_array.append(temp[i])
#            else:
#                temp_array.append(0)
#                
#        temp_array.append(proportions[track])
#        track += 1
#        matrixForML.append(temp_array)
#
#    with open(csvfile, "a") as f:
#        writer = csv.writer(f)
#        writer.writerows(matrixForML)
            
        
def simulation(gene, numOfIter, originalPath, simulation_result_folder, coverage, bt, samtools, art):
    ''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Defining some parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '''
    #Record true variants and their fractions for all simulations
    true_ratios_list = []
    true_variants_list = []
    #Statistics list to output
    totalVarDist_count = []
    totalVarDist_bayes = []
    precision_list = []
    recall_list = []
    pred_object_vals = []
    true_Objective_vals = []
    diff_obj_vals = []
    numOfOptimalSol = list()
    #Some counts
    predictedCorrect_count = 0
    predCorrect_bool_list = []
    seed = 1994
    random.seed(seed)
    
    #Handling some output files
#    outputFolderPath = "{0}/{1}/".format(originalPath, simulation_result_folder)
    outputResultTxtFile = "{0}/{1}/{2}_output_stats.txt".format(originalPath, simulation_result_folder, gene)
    sys.stdout = open(outputResultTxtFile, "w")        #Write print codes to outputResultTxtFile
    
    #Count the number of variants for this gene
    variantsTxtPath = "{0}/sim_data/{1}/variants.txt".format(originalPath,gene)
    num_variants = sum(1 for line in open(variantsTxtPath))
    
    '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Simulation starts here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    
    iteration = 1
    numSimHavingMultSol = 0
    multSolCorrect = 0
    while(iteration < numOfIter+1):  
        true_prop = dict()#dictionary to store the true proportions
        k =random.randint(2,7) #generate a random integer k between 2 and 7
        #generate k random fractions that sum up to 1
        fractions = [random.random() for j in range(k)] 
        s = sum(fractions)
        fractions = [ i/s for i in fractions ]
        #print fractions
        true_variants = []
        #randomly select variants to use for simulated data
        randomVarIndex = random.sample(xrange(1,num_variants+1), k) #start from 1 to num_variants+1 because 0 index of linecache is ''
        for index in randomVarIndex:
            variant =linecache.getline(variantsTxtPath,index) #extract the random variant
            #print variant
            variant = str(variant) 
            variant = variant.rstrip() #remove the "\n" character that is returned by bash
            string1= variant.split(">")
            true_variants.append(string1[1]) #append the variant to the list          
        #print num
        #print true_variants
        
        '''======== Generate the reads using ART ========'''
        total_variants_sequence = ''
        for i in range(len(fractions)):
            variant = true_variants[i]
            simulation_name = true_variants[i] + '_' + str(iteration)+'_'
            file_name = simulation_name + '_reference.fa'
            covOfThisVar = math.ceil((fractions[i]*coverage)) #compute the number of reads to generate
            covOfThisVar = int(covOfThisVar)
            variant_sequence = sh.grep(variant,"{0}/sim_data/{1}/linear.txt".format(originalPath, gene),"-w","-A1") #use bash to extract the variant sequence
            variant_sequence = variant_sequence.rstrip() #remove the "\n" character that is returned by bash
            variant_sequence = str(variant_sequence)
            total_variants_sequence += variant_sequence
            
            #Write sequence to file
            with open(file_name, "w") as sequence_file: 
                sequence_file.write(variant_sequence)
            
            #Set the ART command, I have included a random seed for reproducibility, and a coverage parameter
            ART_command = art+" -qL 33 -qs 10 -qs2 15 -k 3 -rs {} -q -ss HS25 -sam -i ".format(seed) +file_name+" -p -l 76 -f "+str(covOfThisVar)+" -m 200 -s 10 -o "+simulation_name + ' >/dev/null 2>&1'
            os.system(ART_command)
        
        #Putting the pairs together for all variants
        appendFirst_cmd = "cat *_1.fq > "+str(upperfirst(gene)) + "_"+str(iteration)+"_1.fa" #append all the first of the pairs together
        appendSecond_cmd ="cat *_2.fq > "+str(upperfirst(gene))+"_"+str(iteration)+"_2.fa" #append all the second of the pairs together
        os.system(appendFirst_cmd)
        os.system(appendSecond_cmd)
#        ref = upperfirst(gene)+"_bowtie2"
        ref = upperfirst(gene) + "_bowtie"
#        mapping_cmd = 'bowtie2 -x {0} -a -p 4 -1 ./{1}_{2}_1.fa -2 ./{1}_{2}_2.fa -S {1}_{2}.sam >/dev/null 2>&1'.format(ref, upperfirst(gene), str(iteration))
        mapping_cmd = bt+"bowtie -a --best --strata -v 3 -p 4 {0} -1 ./{1}_{2}_1.fa -2 ./{1}_{2}_2.fa --sam {1}_{2}.sam >/dev/null 2>&1".format(ref, upperfirst(gene), str(iteration))

        #Execute commands for bowtie mapping
        os.system(mapping_cmd)
        
        #convert from sam to bam file
        mapped_cmd = samtools+" view -h -F4 {0}_{1}.sam > {0}_{1}_mapped.sam".format(upperfirst(gene),str(iteration))
        paired_cmd = samtools+" view -F8 {0}_{1}_mapped.sam > {0}_{1}_pairedNoHeader.sam".format(upperfirst(gene),str(iteration))
#        singleton_cmd = "samtools view -f8 {0}_{1}_mapped.sam > {0}_{1}_singletonNoHeader.sam".format(upperfirst(gene),str(iteration))
        
        #run the commands
        os.system(mapped_cmd)
        os.system(paired_cmd)
#        os.system(singleton_cmd)
        
        #Tabulate as reads.txt file
        #For paired
        writeReadTable(upperfirst(gene), str(iteration), "paired")
        #For singleton
#        writeReadTable(upperfirst(gene), str(iteration), "singleton")
        
        #Remove unneccessary files for the next iteration.    
        os.system("rm {}*".format(gene))
        os.system("rm *.sam")

        #Keep track of true variants
        true_ratios_list.append(fractions)
        true_variants_list.append(true_variants)
        for j in range(0,len(true_variants)):
            key = true_variants[j]
            true_prop[key] = float(fractions[j])
        
        #Create data matrix where rows=reads and columns=variants
        paired_readsTxt_path = upperfirst(gene)+ '_'+str(iteration)+'_paired_reads.txt'
#        singleton_readsTxt_path = upperfirst(gene)+ '_'+str(iteration)+'_singleton_reads.txt'
        pairedDF = pf.returnMismatchMatrix(paired_readsTxt_path, "paired")
#        singletonDF = returnMismatchMatrix(singleton_readsTxt_path, "singleton")
        dataMatrix = pairedDF
        dataMatrix = dataMatrix.fillna(-1)
        paired_Qmatrix = pf.returnQualityMatrix(paired_readsTxt_path, "paired")
#        singleton_Qmatrix = returnQualityMatrix(singleton_readsTxt_path, "singleton")
        Qmatrix = paired_Qmatrix
        
        #Run the ILP solver
        pred_object_val,var_predicted,reads_cov,all_solutions, all_objective = varSolver.solver(dataMatrix, Qmatrix, "paired")
        #pred_object_val,var_predicted,reads_cov,all_solutions, all_objective = varSolver.solver(dataMatrix)
        #print all_solutions
        #print dataMatrix
#        if len(all_solutions) == 1:
#            continue
        
        
        '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Statistics and Calculations start here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
        #score_list = list()
        #min_score = sys.maxint
        #Qscore_list = [] #list to hold all the average Qsocres for the solutions
        
        #Compute negative log likelihood score for each solution
        #for i in range(len(all_solutions)):
#            print("Solution:{}".format(all_solutions[i]))
#            print("Objective value: {}".format(all_objective[i]))
#            print("Proportions:{}".format(compute_proportions(dataMatrix.loc[reads_cov, all_solutions[i]])))
        #    score = pf.compute_likelihood(dataMatrix.loc[reads_cov, all_solutions[i]], 6)
        #    score_list.append(score)
            
            
        #    if score <= min_score:
        #        min_score = score
                
        #    Qscore_list.append(pf.compute_QSum(Qmatrix.loc[reads_cov,all_solutions[i]]))
            
        #min_Qscore = np.argmin(Qscore_list)
        #var_predicted = all_solutions[min_Qscore]
        #minQIndices = np.where(np.array(Qscore_list) == np.array(Qscore_list).min())[0]
        
        #if len(minQIndices) > 1:
        #    print("More than 1 solution having minimum quality score")
            
        #Keep track of precision and recall for all simulations
        precision_list.append(precision(var_predicted, true_variants))
        recall_list.append(recall(var_predicted, true_variants))
        pred_object_vals.append(pred_object_val)
        
        #Construct dataframe of true variants and predicted variants
        #true_DF = dataMatrix.loc[reads_cov,true_variants]
#        true_DF =  true_DF[(true_DF.T != -1).any()]
        #predicted_DF = dataMatrix.loc[reads_cov,var_predicted]
        predicted_DF = Qmatrix.loc[reads_cov,var_predicted]
        prop_count = pf.compute_proportions(predicted_DF)
        #prop_bayes = bayes_compute_proportions(predicted_DF)
#        if proportion_method == "count":
#            prop = count_compute_proportions(predicted_DF)
#        elif proportion_method == "bayes":
#            prop = bayes_compute_proportions(predicted_DF)
        pred_prop_count = pf.create_dictionary(var_predicted, prop_count)
        #pred_prop_bayes = create_dictionary(var_predicted, prop_bayes)
        val_count = 100.0*totalVariationDist(pred_prop_count, true_prop)
        #val_bayes = totalVariationDist(pred_prop_bayes, true_prop)
        totalVarDist_count.append(val_count)
        #totalVarDist_bayes.append(val_bayes)
        
        #Print true and predicted variants
        print "======================================== SIMULATION " + str(iteration) + " ====================================================" +  "\n"
        print("Number of optimal solutions: {}".format(len(all_solutions)))
        numOfOptimalSol.append(len(all_solutions))
        print all_solutions
        #print ('Qscores are: {}').format(Qscore_list)
        #print("Minimum q score is :{}".format(Qscore_list[min_Qscore]))
#        print ('Likelihood scores are: {}').format(score_list)
#        print("Minimum likelihood score is :{}".format(score_list[minLikeli]))
        print("True variants are: {}\n".format(true_variants))
        print("Predicted variants are: {}\n".format(var_predicted))
        print("True proportions are: {}\n".format(true_prop))
        #print("Predicted proportions using Bayes method are: {}\n".format(pred_prop_bayes))
        print("Predicted proportions using Counting method are: {}\n".format(pred_prop_count))
        #correct = False
        #if set(true_variants) == set(var_predicted):
        #    correct = True
            
        #if correct:
        #    print("This sim is predicted correctly based on quality score.")
        #else:
        #    print("This sim is not predicted correctly based on quality score.")
        
        if len(all_solutions) > 1:
            numSimHavingMultSol += 1
            
            #if correct:
            #    multSolCorrect += 1
        
        #Compute objective value of true variants
        #true_Objective_val, bad_reads = compute_true_objective_val(true_DF)
        #true_Objective_vals.append(true_Objective_val)
        
        #Compute the difference in objective vlaue: Predicted - True
        #diff_obj_vals.append(compute_obj_diff(pred_object_val,true_Objective_val))
        
        #Count simulations in which the ILP predicted correctly
        if predictedCorrectly(var_predicted, true_variants):
                predictedCorrect_count += 1
                predCorrect_bool_list.append(True)
        else:
                predCorrect_bool_list.append(False)
                
        #if len(bad_reads) != 0:
        #    print(dataMatrix.loc[bad_reads,:])   
            
        iteration += 1

    
    '''~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Here is the summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'''
    
    print "======================================== {0}: SUMMARY STATISTICS ====================================================\n".format(gene)
#    avg_totalVarDist_bayes = sum(totalVarDist_bayes)/len(totalVarDist_bayes)
#    true_avg_totalVarDist_bayes = sum(list(itertools.compress(totalVarDist_bayes, predCorrect_bool_list)))/sum(predCorrect_bool_list)
#    variance_totalVarDist_bayes = map(lambda x: (x - avg_totalVarDist_bayes)**2, totalVarDist_bayes)
#    true_variance_totalVarDist_bayes = map(lambda x:(x - true_avg_totalVarDist_bayes)**2, list(itertools.compress(totalVarDist_bayes, predCorrect_bool_list)))
#    variance_totalVarDist_bayes = sum(variance_totalVarDist_bayes)/len(variance_totalVarDist_bayes)
#    true_variance_totalVarDist_bayes = sum(true_variance_totalVarDist_bayes)/len(true_variance_totalVarDist_bayes)
#    std_totalVarDist_bayes = math.sqrt(variance_totalVarDist_bayes)
#    true_std_totalVarDist_bayes = math.sqrt(true_variance_totalVarDist_bayes)
    
    avg_totalVarDist_count = sum(totalVarDist_count)/len(totalVarDist_count)
    true_avg_totalVarDist_count = sum(list(itertools.compress(totalVarDist_count, predCorrect_bool_list)))/sum(predCorrect_bool_list)
    variance_totalVarDist_count = map(lambda x: (x - avg_totalVarDist_count)**2, totalVarDist_count)
    true_variance_totalVarDist_count = map(lambda x:(x - true_avg_totalVarDist_count)**2, list(itertools.compress(totalVarDist_count, predCorrect_bool_list)))
    variance_totalVarDist_count = sum(variance_totalVarDist_count)/len(variance_totalVarDist_count)
    true_variance_totalVarDist_count = sum(true_variance_totalVarDist_count)/len(true_variance_totalVarDist_count)
    std_totalVarDist_count = math.sqrt(variance_totalVarDist_count)
    true_std_totalVarDist_count = math.sqrt(true_variance_totalVarDist_count)
    
    context=1
    print "({0})Total Variation Distance:\n".format(context)
    print "Counting Method ~ Total variation distances are:",totalVarDist_count, "\n"
    #print "Bayes' Method ~ Total variation distances are:",totalVarDist_bayes, "\n"
    print "Counting Method ~ The average of total variation distance is:", avg_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The average of total variation distance is:", avg_totalVarDist_bayes, "\n"
    #print "The variance of total variation distance is:", variance_totalVarDist, "\n"
    print "Counting Method ~ The standard deviation of total variation distance is:",std_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The standard deviation of total variation distance is:",std_totalVarDist_bayes, "\n"
    context+=1
    
    print "({0})Total Variation Distance for variants which are predicted correctly:\n".format(context)
    print "Counting Method ~ Total variation distances are:",list(itertools.compress(totalVarDist_count, predCorrect_bool_list)), "\n"
    #print "Bayes' Method ~ Total variation distances are:",list(itertools.compress(totalVarDist_bayes, predCorrect_bool_list)), "\n"
    print "Counting Method ~ The average of total variation distance is:", true_avg_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The average of total variation distance is:", true_avg_totalVarDist_bayes, "\n"
    #print "The variance of total variation distance is:", variance_totalVarDist, "\n"
    print "Counting Method ~ The standard deviation of total variation distance is:",true_std_totalVarDist_count, "\n"
    #print "Bayes' Method ~ The standard deviation of total variation distance is:",true_std_totalVarDist_bayes, "\n"
    context+=1
    
    avg_prec = sum(precision_list)/len(precision_list)
    std_prec = np.std(np.array(precision_list))
    print "({0}) Precision: \n".format(context)
    print 'Precision is:', precision_list, "\n"
    print "Average of precision is: ", avg_prec, "\n"
    print "Standard deviation of precision is: ", std_prec, "\n"
    context+=1
    
    avg_rec = sum(recall_list)/len(recall_list)
    std_rec = np.std(np.array(recall_list))
    print "({0}) Recall : \n".format(context)
    print 'Recall is:', recall_list, "\n"
    print "Average of recall is: ", avg_rec, "\n"
    print "Standard deviation of recall is: ", std_rec, "\n"
    context+=1
    
    #avg_diffObjVal = sum(diff_obj_vals)/len(diff_obj_vals)
    #std_diffObjVal = np.std(np.array(diff_obj_vals))
    #print "({0}) Objective Value: \n".format(context)
    #print 'Predicted objective values are:', pred_object_vals, "\n"
    #print 'True objective values are:', true_Objective_vals, "\n"
    #print 'The difference in objective values are:', diff_obj_vals, "\n"
    #print "Average of difference in objective value is: ", avg_diffObjVal, "\n"
    #print "Standard deviation of difference in objective value is: ", std_diffObjVal, "\n"
    #context+=1
    
    print "({0})Accuracy: \n".format(context)
    print 'Total number of simulations: ', numOfIter , "\n"
    print("Number of simulations which are predicted correctly: {0}\n".format(predictedCorrect_count))
    print 'Percentage of simulations predicted correctly: ', 100*predictedCorrect_count/numOfIter, "%\n"
    context+=1
    
    print "({0})Optimal solutions: \n".format(context)
    print("Number of optimal solutions: {0}\n".format(numOfOptimalSol))
    print("Average number of optimal solutions: {0}\n".format(np.mean(numOfOptimalSol)))
    context+=1
    
    #print "({0})Statistics about quality score: \n".format(context)
    #print("Number of simulations having multiple solutions: {}".format(numSimHavingMultSol))
    #print("Number of these simulations which are correct: {}".format(multSolCorrect))
    #print("The percentage is : {} %".format(100.0*(multSolCorrect/numSimHavingMultSol)))
    #context+=1
    
    sys.stdout.close()
    #return precision_list, recall_list, diff_obj_vals, totalVarDist_count
    return precision_list, recall_list, totalVarDist_count
