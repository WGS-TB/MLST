#!/usr/bin/python
from __future__ import division
from collections import defaultdict
from scipy.special import comb
import pandas as pd
import sh
import math
import numpy as np
import random
import os
import argparse
import matplotlib.pyplot as plt
import variantILP as ilp
import itertools
import sys
import cplex

NO_BINOM = False
TEST_EMPTY_LIST = True

# pred_prop and true_prop: dictionaries
def totalVariationDist(pred_prop, true_prop):
    common_keys = set.intersection(set(pred_prop.keys()), set(true_prop.keys()))
    diff_keys_pred = set.difference(set(pred_prop.keys()), common_keys)
    diff_keys_true = set.difference(set(true_prop.keys()), common_keys)
    common_keys = list(common_keys)
    diff_keys_pred = list(diff_keys_pred)
    diff_keys_true = list(diff_keys_true)
    
    totalVarDist=0
    for key in common_keys:
        totalVarDist += abs(pred_prop[key] - true_prop[key])
        
    for key in diff_keys_pred:
        totalVarDist += pred_prop[key]
        
    for key in diff_keys_true:
        totalVarDist += true_prop[key]
        
    return totalVarDist/2

def upperfirst(x):
    return x[0].upper() + x[1:]

def precision(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(predicted))
        
def recall(predicted, true):
    truePos = set.intersection(set(predicted), set(true))
    return float(len(truePos)/len(true))

def tree():
    return defaultdict(tree)

def predictedCorrectly(predicted, true):
        return set(predicted) == set(true)

def solver(dataMatrix):
#    data_matrix = returnDataMatrix("/home/glgan/Documents/Borrelia/data/simData/clpA_7_weighted.csv")
#    data_matrix = returnDataMatrix(dataFile+fileName)
    data_matrix = dataMatrix
    total_read = data_matrix.shape[0]
    total_var = data_matrix.shape[1]  #As first column are reads
    
    xIsVariant = pd.DataFrame(index=[var for var in data_matrix.columns.tolist()], data=["x{0}".format(i+1) for i in range(total_var)], columns=["Variable"])
    variantName_variable_dict = xIsVariant["Variable"].to_dict()
#    variable_variantName_dict = pd.DataFrame(data=[var for var in data_matrix.columns.tolist()[1:]], index=["x{0}".format(i+1) for i in range(total_var)], columns=["Variant"])["Variant"].to_dict()
    reads = data_matrix.index.tolist()
    readAndMM = list()
    
    #Grab unique set of mismatches for each read 
    for row in data_matrix.itertuples():
        varMM = list(set(list(row[1:])))
        varMM = [int(mm) for mm in varMM if mm != -1]
        readAndMM.append(varMM)
            
    readAndMMVariable = [[[j,"y{0}_{1}".format(i+1, j)] for j in readAndMM[i]] for i in range(total_read)]
    yIsRead = pd.DataFrame(data=reads, columns=["Read"])
    yIsRead["Variable"] = readAndMMVariable
    
    #Dictionary for easier retrievable of y_variable name when formulating ILP
    #key is a tuple of read name and number of mismatches, value = read variable name
    readMM_varName_dict = dict()
    #key is y variable name, value is the matching read
    varName_read_dict = dict()
    for r in range(total_read):
        for mm in readAndMM[r]:
            readMM_varName_dict[(reads[r], mm)] = "y{0}_{1}".format(r+1, mm)
            varName_read_dict["y{0}_{1}".format(r+1, mm)] = reads[r]
            
           
    '''=============================================== Form ILP ====================================================== '''
    model = cplex.Cplex()
    model.objective.set_sense(model.objective.sense.minimize)
    #add variables related to variants, x_j represents variant j
    model.variables.add(obj=[1]*total_var, names=xIsVariant["Variable"].tolist(), types=[model.variables.type.binary]*total_var)
    #model.variables.add(names=xIsVariant["Variable"].tolist(), types=[model.variables.type.binary]*total_var)
    
    #add variables related to reads, y_ik means read i with k mismatches
    y_weights = [mm[0] for mmName in readAndMMVariable for mm in mmName]
    y_variables = [mm[1] for mmName in readAndMMVariable for mm in mmName]
    model.variables.add(obj=y_weights, names=y_variables, types=[model.variables.type.binary]*len(y_variables))
    
    #Constraint: Cover most of the read
    keep_alpha = 1.0
    model.linear_constraints.add(lin_expr=[[y_variables, [1]*len(y_variables)]], senses=["G"], rhs=[keep_alpha*total_read], names=["c{0}".format(1+model.linear_constraints.get_num())])
    
    #Constraint: For each read i, only allow it to be covered with a unique number of mismatch i.e. 
    #only one of y_i0, y_i1, y_i2, y_i3 is =1
    uniqueMMConstr = list()
    for oneRead in readAndMMVariable:
        temp = list()
        yVariable = list()
        for i in range(len(oneRead)):
            yVariable.append( (oneRead[i])[1] )
        temp.append(yVariable)
        temp.append([1]*len(yVariable))
        uniqueMMConstr.append(temp)
        
    model.linear_constraints.add(lin_expr=uniqueMMConstr, senses=["L"]*total_read, rhs=[1]*total_read, names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(total_read)])
    
    #Limit maximum number of variants
    #model.linear_constraints.add(lin_expr=[[xIsVariant["Variable"].tolist(), [1]*total_var]], senses=["L"], rhs=[10], names=["c{0}".format(1+model.linear_constraints.get_num())])
    
    #This is for easier formulation of the constraints related to: If y_ik=1, then it must be 
    #covered by some variant with k mm
    yCoverConstr = list()
    readVar_variant_dict = dict()   #key is read variable name, value is the list of variants covering that read(with the particular mm)
    for row in data_matrix.itertuples():
        read = row[0]
        uniqueMm = list(set(row[1:]))
        uniqueMm = [int(mm) for mm in uniqueMm if mm!=-1]
        
        for mm in uniqueMm:
            temp = [i for i in range(len(list(row[1:]))) if list(row[1:])[i] == mm]
            
            xVariable = [variantName_variable_dict[variant] for variant in data_matrix.columns[temp].tolist()]
            yVariable = [ readMM_varName_dict[read, mm] ]
            readVar_variant_dict[readMM_varName_dict[read, mm]] = xVariable
            yCoverConstr.append([ xVariable + yVariable, [1]*len(xVariable) + [-1]*len(yVariable) ])
        
    #Constraint: If y_ik is chosen, it must be covered by some variant j with k mm
    model.linear_constraints.add(lin_expr=yCoverConstr, senses=["G"]*len(y_variables), rhs=[0]*len(y_variables), names=["c{0}".format(i+1+model.linear_constraints.get_num()) for i in range(len(y_variables))])
    
    #model.write("a.lp")
    model.set_results_stream(None)
#    model.solve()
    
    #options for searching more optimal solutions
    #model.parameters.mip.pool.capacity.set(10)
    model.parameters.mip.pool.intensity.set(4)
    #model.parameters.mip.limits.populate.set(2100000000)
    model.parameters.mip.pool.absgap.set(0)
    model.parameters.mip.pool.replace.set(1)
    model.populate_solution_pool()
       
    '''=========================================== Presenting Results ========================================================'''
    objvalue = model.solution.get_objective_value()
    varNames = model.variables.get_names()
    varValues = model.solution.get_values(varNames)
    conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
    conclusion["Decision Variable"] = varNames
    conclusion["Value"] = varValues

    present = conclusion[conclusion["Value"]==1]
    variantsPresent = xIsVariant[xIsVariant["Variable"] .isin(present["Decision Variable"].tolist())]
    varPresentList = variantsPresent.index.tolist()
#    xPresentList = variantsPresent["Variable"].tolist()

    yCovered = present[present["Decision Variable"].str.contains(u"y.*")]["Decision Variable"].tolist()
    readsCovered = [varName_read_dict[y] for y in yCovered]
    
    allSol = list()
    allObjValue = list()
    for i in range(model.solution.pool.get_num()):
        objvalue = model.solution.pool.get_objective_value(i)
        varNames = model.variables.get_names()
        varValues = model.solution.pool.get_values(i,varNames)
        conclusion = pd.DataFrame(columns=["Decision Variable", "Value"])
        conclusion["Decision Variable"] = varNames
        conclusion["Value"] = varValues
    
        present = conclusion[conclusion["Value"]==1]
        variantsPresent = xIsVariant[xIsVariant["Variable"] .isin(present["Decision Variable"].tolist())]
        varPresentList = variantsPresent.index.tolist()
        allSol.append(varPresentList)
        allObjValue.append(objvalue)
    
    return objvalue, varPresentList, readsCovered, allSol, allObjValue

#function to generate matrix to be passed to 
def Generate_Matrix(path):
    var_list = [] #holds the variants
    read_list = [] #holds the reads
    mismatch_list = [] #holds the mismatches
    with open(path) as inf:
        for line in inf:
            parts = line.split('\t') # split line into parts
            if len(parts) > 1:   # if at least 2 parts/columns
                read_list.append(parts[0]) #append reads to a list
                var_list.append(parts[1])   #append vars to a list
                mismatch_list.append(parts[2]) #append mismatches to list
        flag = True #makes sure all the previous steps completed successfully
    if flag is True:
        read_var_dict = defaultdict(list) #dictionary holding all the variants that a read maps to
        read_mismatch_dict = defaultdict(list) #dictionary holding all the mismatches that a read has to its variants
        d_2 = defaultdict(list) #dictionary holding indices for later use


        for i in range(len(read_list)):
            num_mismatch = mismatch_list[i].count('>') #count the number of mismatches for each read
            #append the appropraite suffix for paired reads
            '''if  i%2 == 0:    
                read_list[i] = read_list[i]+'-1/2'
            else:
                read_list[i] = read_list[i]+'-2/2'
            '''
            read_var_dict[read_list[i]].append(var_list[i]) #append all the variants that read read_i maps to
            read_mismatch_dict[read_list[i]].append(num_mismatch) #append all the mismatches that each read_i has when it maps to a variants
            d_2[read_list[i]].append(i) #for testing purposes

        var_list = set(var_list) #removes duplicates
        matrix_dict = tree() #creates a 2-D dictionary object later used to generate the read-variant matrix datastructure
 
#       create a 2D dictionary that contains all the possible combinations of a read with a variant and the number of mismatches.
        for var in var_list:
            for read in read_var_dict:   #key=read name
                temp_var_list = read_var_dict[read]   #list of variants that key maps to 
                if var in temp_var_list:
                    index = temp_var_list.index(var)
                    mismatch = read_mismatch_dict[read][index] #get the number of mismatches
                    #only for perfect matches.
                    #if val == 1:
                    #print val
                    matrix_dict[read][var] = int(mismatch) #add it to the matrix data structure

        df = pd.DataFrame(matrix_dict).T.fillna(-1) #convert 2-D dictionary to a matrix
    return df

def compute_probability(n, k):
    #compute the probability of read of length n mapping to a variant with k mismatches using the binomial distribution
    if NO_BINOM:
        b=10**10
    else:
        b = comb(n, k, exact=False)
    
    x = math.pow(0.99,(n-k))
    y = math.pow(0.01,k)
    prob = b*x*y
        
    return prob

def compute_proportions(dataframe):
    #computes the proportion of a set of variants given a set of reads uing probabilistic methods
        prob_list = [] #a list to hold the probabilities
        for row in dataframe.itertuples(index=False):
                temp_list = list(row)
                #compute the probability for each row in the matrix
                for i in range(len(temp_list)):
                    if temp_list[i] >= 0:
                        temp_list[i] = compute_probability(76,int(temp_list[i]))
                    else:
                        temp_list[i] = 0
                total = sum(temp_list)
                #solve for k
                try:
                    temp_list = [j*(1.0/total) for j in temp_list]
                except ZeroDivisionError:
                    print(total)
                    print(temp_list)
                prob_list.append(temp_list)
        col_sums = [sum(k) for k in zip(*prob_list)]
        total_sum = sum(col_sums)
        prop_list = [100.0*l*(1/total_sum) for l in col_sums]
        return prop_list     

def compute_True_objective_val(dataframe):
        run_sum = 0
        bad_read = list()
        #Even no limit on mm, there are reads which are covered by predicted variants but not
        #true variants. Identify largest mm , if this case happens, increment by max_mm+1
        max_mm = max(dataframe.max())
        for row in dataframe.itertuples():
                my_list = [i for i in list(row)[1:] if i >= 0]
                if TEST_EMPTY_LIST == True:
                    if len(my_list) > 0:
                        run_sum += min(my_list)
                    else:
                        run_sum+=max_mm + 1
                        bad_read.append(list(row)[0])
                        print(list(row))
                else:
                    if len(my_list) == 0:
                        bad_read.append(list(row)[0])
                        print(list(row))
                    run_sum += min(my_list)
        run_sum += len(dataframe.columns)
        return run_sum, bad_read

def create_dictionary(keys, vals):
        my_dict = dict()
        if len(keys) == len(vals):
                for i in range(len(keys)):
                        my_dict[keys[i]] = vals[i]
        return my_dict 

def Compute_obj_diff(predicted, true):
        diff = predicted - true
        return diff
    
def compute_likelihood(df):
    numVar = df.shape[1]
    likelihood_list = list()
    max_mm = 2
    
    for row in df.itertuples(index=False):
        read = list(row)
        
        temp = list()
        for i in range(numVar):
            if read[i] == -1:
                prob = (0.01)**(max_mm+1) * (0.99)**(76 - max_mm -1)
                temp.append(prob)
            else:
                prob = (0.01)**(read[i]) * (0.99)**(76 - read[i])
                temp.append(prob)
                
        likelihood_list.append( sum(temp) )
    
    likelihood_list = [i/(2.0*76*numVar) for i in likelihood_list]
    neg_log_likelihood = [-1.0*np.log10(j) for j in likelihood_list]
    
    score = sum(neg_log_likelihood)
    return score

ap = argparse.ArgumentParser()
ap.add_argument("-g", "--gene", required = True,  help="name of gene")
ap.add_argument("-l", "--numOfIter", required = False, default = 40, type=int)
ap.add_argument("-o", "--outputFolder", required=True)
args = vars(ap.parse_args())
random.seed(a=1994) #set random seed
true_ratios = []
true_variants = []
total_var = []
Precision = []
Recall = []
pred_object_vals = []
true_Objective_vals = []
diff_obj_vals = []
count = 0
bool_list = []
minNegLogLike_correct = 0
minSizeOpt_count = 0
minNegLogLike_hasTrue = 0

for x in range(1,args["numOfIter"]+1):  
    k =random.randint(2,7) #generate a random integer k between 2 and 7
    #generate k random fractions that sum up to 1
    fractions = [random.random() for j in range(k)] 
    s = sum(fractions)
    fractions = [ i/s for i in fractions ]
	#run a bash sub command to sort the variants text file and return the first k variants
    variants = (sh.head(sh.sort("variants.txt", "-R"), "-n", k))
    variants = list(variants) #convert output from runningCommand to a list
    #print variants #for testing purposes
    ratios = [] #list to store the proprotions of reads to generate
    true_prop = dict()#dictionary to store the true proportions
    gene = args["gene"] #get the locus/gene of interest
    gene = upperfirst(gene)
    total_variant_seqs = "" #string for all the current variant sequences
    file_name2 = gene+"_"+str(x)+".fas"
    variants_current = [] #list to hold the variants randomly chosen
    for i in range(len(fractions)):        
        variant = str(variants[i]) 
        variant = variant.rstrip() #remove the "\n" character that is returned by bash
        string1= variant.split(">")
        variants_current.append(string1[1]) #append the variant to the list
        sim_name = string1[1]+"_"+str(x)+"_" #name of simulation
        file_name = sim_name+"reference.fa"
        num_reads = math.ceil(fractions[i]*197) #compute the number of reads to generate for each variant
        ratios.append(int(num_reads)) 
        variant_seq = sh.grep(variant,"linear.txt","-w","-A1") #use bash to extract the variant sequence
        variant_seq = variant_seq.rstrip() #remove the "\n" character that is returned by bash
        variant_seq = str(variant_seq)
        total_variant_seqs = total_variant_seqs + variant_seq
        #print temp2 
        cmd = "art_illumina -q -ss HS25 -sam -i "+file_name+" -p -l 76 -c "+str(ratios[i])+" -m 200 -s 10 -o "+sim_name + " >/dev/null"#the ART command to generate the simulated data for a variant
        #write the variant sequence to a text file
        with open(file_name, "w") as text_file: 
            text_file.write(variant_seq)
                #print sim_name #testing purposes
                #print file_name #testing purposes
                #print cmd #testing purposes
        os.system(cmd) #run the ART command
    with open(file_name2, "w") as text_file2:
		text_file2.write(total_variant_seqs)
    new_cmd = "cat *_1.fq > "+gene+"_"+str(x)+"_1.fa" #append all the first of the pairs together
    new_cmd2 ="cat *_2.fq > "+gene+"_"+str(x)+"_2.fa" #append all the second of the pairs together
    ref = upperfirst(gene)+"_bowtie"
    new_cmd3 = "bash /home/glgan/Documents/Borrelia/scripts/temp.sh "+ gene+"_"+str(x) + " " + gene + "_" + str(x)+ " " + ref + " >/dev/null "
    #print new_cmd3 
    os.system(new_cmd) #run the command
    os.system(new_cmd2) #run the command
    os.system(new_cmd3)
    os.system("rm "+args["gene"]+"*") #remove unneccessary files for the next iteration.
    #print ratios #for testing purposes
    X =  gene+"_"+str(x)
    X = upperfirst(X)
    #print X
    
    true_ratios.append(fractions)
    true_variants.append(variants_current)
    for j in range(0,len(variants_current)):
            key = variants_current[j]
            true_prop[key] = float(fractions[j])*100
    
    path = X+'_reads.txt'
    df = Generate_Matrix(path)
    df.rename(columns={'Unnamed: 0': 'Read'}, inplace=True)
    pred_object_val,var_predicted,reads_cov,all_solutions, all_objective = solver(df)
    
#    If there is one solution among all optimal which matches true variants, assign to var_predicted
    for solution in all_solutions:
        if set(solution) == set(variants_current):
            var_predicted = solution
            break
    
    print "======================================== SIMULATION " + str(x) + " ====================================================" +  "\n"    
    #Likelihood approach
    score_list = list()
    min_score = sys.maxint
    print("Number of optimal solutions: {}".format(len(all_solutions)))
    
    #compute negative log likelihood score for each solution
#    for i in range(len(all_solutions)):
#        print("Solution:{}".format(all_solutions[i]))
#        print("Objective value: {}".format(all_objective[i]))
#        print("Proportions:{}".format(compute_proportions(df.loc[reads_cov, all_solutions[i]])))
#        score = compute_likelihood(df.loc[reads_cov, all_solutions[i]])
#        score_list.append(score)
#        
#        if score <= min_score:
#            min_score = score
        
#        print("\nNegative log likelihood score:{}\n".format(score))
    
    #identify those solutions which have minimum negative log likelihood
#    min_sol_list = [all_solutions[i] for i in range(len(all_solutions)) if score_list[i] == min_score]
#    var_predicted = min_sol_list[0]
#    similarity=0.001
#    similar_likelihood_sol_list = [all_solutions[i] for i in range(len(all_solutions)) if (score_list[i] <= (1+similarity)*min_score)]
#    similar_likelihood_score_list = [score_list[i] for i in range(len(all_solutions)) if (score_list[i] <= (1+similarity)*min_score)]
#    var_predicted = list(set().union(*similar_likelihood_sol_list))
#    print("Similar likelihood list: {}".format(similar_likelihood_sol_list))
#    print("Similar likelihood score list: {}".format(similar_likelihood_score_list))
    
    Precision.append(precision(var_predicted, variants_current))
    Recall.append(recall(var_predicted, variants_current))
    pred_object_vals.append(pred_object_val)
    df1 = df[var_predicted]
    df2 = df.loc[reads_cov,variants_current]
    df3 = df.loc[reads_cov,var_predicted]
    prop = compute_proportions(df3)
    pred_prop = create_dictionary(var_predicted, prop)
    val = totalVariationDist(pred_prop, true_prop)
    total_var.append(val)

    print 'True variants are:', variants_current, "\n"
    print 'Predicted variants are:', var_predicted, "\n"
    print 'True proportions are:', true_prop, "\n"
    print 'Predicted proportions are:', pred_prop, "\n"
    #print 'Solution(s):', all_solutions, "\n"
    
    #Observe whether minimum negative log likelihood really gives the true variants
#    for sol in min_sol_list:
#        if set(sol) == set(variants_current):
#            minNegLogLike_correct += 1
#            min_sol = sol
#            break
    
    #Observe whether minimum negative log likelihood solutions implies minimum number of variants
#    minVariants = min(map(len,all_solutions))
#    if len(min_sol) == minVariants:
#        minSizeOpt_count +=1
        
    #Count number of simulations where solution with minimum negative log likelihood contains true variants
#    for sol in min_sol_list:
#        if len(set(variants_current) - set(sol)) == 0:
#            minNegLogLike_hasTrue += 1
#            break
        
#    print("Negative log likelihood score list:{}".format(score_list))
#    print("The solution which has minimum negative log likelihood is {0} and its score:{1}".format(min_sol, min_score))
    
    true_Objective_val, bad_reads = compute_True_objective_val(df2)
    true_Objective_vals.append(true_Objective_val)
    diff_obj_vals.append(Compute_obj_diff(pred_object_val,true_Objective_val))
    
    if predictedCorrectly(var_predicted, variants_current):
            count += 1
            bool_list.append(True)
    else:
            bool_list.append(False)
            
    if len(bad_reads) != 0:
        print(df.loc[bad_reads,:])

print "======================================== {0}: SUMMARY STATISTICS ====================================================\n".format(gene)
avg_totalVarDist = sum(total_var)/len(total_var)
true_avg_totalVarDist = sum(list(itertools.compress(total_var, bool_list)))/sum(bool_list)
variance_totalVarDist = map(lambda x: (x - avg_totalVarDist)**2, total_var)
true_variance_totalVarDist = map(lambda x:(x - true_avg_totalVarDist)**2, list(itertools.compress(total_var, bool_list)))
variance_totalVarDist = sum(variance_totalVarDist)/len(variance_totalVarDist)
true_variance_totalVarDist = sum(true_variance_totalVarDist)/len(true_variance_totalVarDist)
std_totalVarDist = math.sqrt(variance_totalVarDist)
true_std_totalVarDist = math.sqrt(true_variance_totalVarDist)
context=1
print "({0})Total Variation Distance:\n".format(context)
print "Total variation distances are:",total_var, "\n"
print "The average of total variation distance is:", avg_totalVarDist, "\n"
#print "The variance of total variation distance is:", variance_totalVarDist, "\n"
print "The standard deviation of total variation distance is:",std_totalVarDist, "\n"
context+=1
print "({0})Total Variation Distance for variants which are predicted correctly:\n".format(context)
print "Total variation distances are:",list(itertools.compress(total_var, bool_list)), "\n"
print "The average of total variation distance is:", true_avg_totalVarDist, "\n"
#print "The variance of total variation distance is:", variance_totalVarDist, "\n"
print "The standard deviation of total variation distance is:",true_std_totalVarDist, "\n"
context+=1

avg_prec = sum(Precision)/len(Precision)
std_prec = np.std(np.array(Precision))
print "({0}) Precision: \n".format(context)
print 'Precision is:', Precision, "\n"
print "Average of precision is: ", avg_prec, "\n"
print "Standard deviation of precision is: ", std_prec, "\n"
context+=1

avg_rec = sum(Recall)/len(Recall)
std_rec = np.std(np.array(Recall))
print "({0}) Recall : \n".format(context)
print 'Recall is:', Recall, "\n"
print "Average of recall is: ", avg_rec, "\n"
print "Standard deviation of recall is: ", std_rec, "\n"
context+=1

avg_diffObjVal = sum(diff_obj_vals)/len(diff_obj_vals)
std_diffObjVal = np.std(np.array(diff_obj_vals))
print "({0}) Objective Value: \n".format(context)
print 'Predicted objective values are:', pred_object_vals, "\n"
print 'True objective values are:', true_Objective_vals, "\n"
print 'The difference in objective values are:', diff_obj_vals, "\n"
print "Average of difference in objective value is: ", avg_diffObjVal, "\n"
print "Standard deviation of difference in objective value is: ", std_diffObjVal, "\n"
context+=1

print "({0})Accuracy: \n".format(context)
print 'Total number of simulations: ', args["numOfIter"] , "\n"
print("Number of simulations which are predicted correctly: {0}\n".format(count))
print 'Percentage of simulations predicted correctly: ', 100*count/x, "%\n"
context+=1

#print "({0})Numbers related to likelihood approach: \n".format(context)
#print 'Number of simulations where solution with minimum negative log likelihood is the true solution: ', minNegLogLike_correct, "\n"
#print('Out of the {0} simulations which are predicted correctly, {1} simulations have optimal solutions which do not have a minimum negative likelihood score.\n'.format(count, count-minNegLogLike_correct) )
#print('The percentage of minimum negative log likelihood score giving true solution: {}\n'.format(minNegLogLike_correct*100.0/count) )
#print 'Number of simulations where solution with minimum negative log likelihood has minimum number of variants: ', minSizeOpt_count, "\n"
#print("Number of simulations where solution with minimum negative log likelihood CONTAINS true variants: {0}\n".format(minNegLogLike_hasTrue))
X = []
os.system("rm {0}_*.fa".format(gene))

for i in range(1,args["numOfIter"]):
        X.append(i)
plt.figure()
plt.hist(total_var, bins=args["numOfIter"]-10)
plt.xlabel('Total Variation Distance in %')
plt.ylabel('Frequency')
#plt.show()
plt.savefig("/home/glgan/Documents/Borrelia/{0}/{1}_totalVarDist".format(args["outputFolder"], gene))

plt.figure()
plt.hist(Precision, bins=args["numOfIter"]-10)
plt.xlabel('Precision')
plt.ylabel('Frequency')
plt.savefig("/home/glgan/Documents/Borrelia/{0}/{1}_precision".format(args["outputFolder"], gene))

plt.figure()
plt.hist(Recall, bins=args["numOfIter"]-10)
plt.xlabel('Recall')
plt.ylabel('Frequency')
plt.savefig("/home/glgan/Documents/Borrelia/{0}/{1}_recall".format(args["outputFolder"], gene))

plt.figure()
plt.hist(diff_obj_vals, bins=args["numOfIter"]-10)
plt.xlabel('Difference in objective values: Predicted - True')
plt.ylabel('Frequency')
plt.savefig("/home/glgan/Documents/Borrelia/{0}/{1}_diffObjVals".format(args["outputFolder"], gene))