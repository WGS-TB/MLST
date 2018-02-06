#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 13:33:45 2017

@author: glgan
"""

import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.spatial.distance import hamming
import itertools

matplotlib.style.use('ggplot')

def generateNewRefNewStr(refStrains, loci, writePath, df_justStrains_new):
    seed=13
    numLoci=len(loci)
    writePath = os.path.abspath(writePath) 
    reference = pd.read_csv(refStrains,sep="\t",usecols=range(1,numLoci+1))

    for i in range(1,11):
        adding = df_justStrains_new.sample(frac=0.5, random_state=seed+i)
        
        for l in loci:
            adding.loc[:,l] = adding.loc[:,l].str.split("_").str[1]
        
        adding.reset_index(drop=True, inplace=True)
        temp = reference.append(adding[loci])
        temp = temp.reset_index(drop=True)
        adding[loci].to_csv("{0}/added_to_ref{1}.csv".format(writePath,i))
        temp.to_csv(os.path.join(writePath, "strain_newRef{}.txt".format(i)), sep="\t")

def generateNewRefRemoveExist(refStrains, remove, loci, ind, writePath):
    numLoci=len(loci)
    reference = pd.read_csv(refStrains,sep="\t",usecols=range(1,numLoci+1))
    writePath = os.path.abspath(writePath) 
    for name in loci:
        reference["%s" %name] = name + "_" + reference["%s" %name].astype(str)
    
    new_ref = reference[~reference[loci].isin(remove).all(axis=1)]
    
    for l in loci:
        new_ref.loc[:,l] = new_ref.loc[:,l].str.split("_").str[1]
    
    new_ref.reset_index(drop=True, inplace=True)
    new_ref.to_csv(os.path.join(writePath, "exist_strain_newRef{}.txt".format(ind)), sep='\t')

def plotDistribution(df_justStrains, loci, name):
    df_onlyIndex = df_justStrains[:]

    for l in loci:
        df_onlyIndex[l] = df_justStrains[l].str.split('_').str[1]
    
    matrix_onlyIndex = df_onlyIndex[loci].as_matrix()
    ST_alleles_dict = {i:ST for i,ST in itertools.izip(range(len(matrix_onlyIndex)), df_onlyIndex.ST.tolist())}
    hammingMatrix = np.zeros( (matrix_onlyIndex.shape[0], matrix_onlyIndex.shape[0]) )
    hammingMatrix_asdict = {(i,j):0.0 for (i,j) in itertools.product(df_onlyIndex.ST.tolist(), df_onlyIndex.ST.tolist())}
    hamming_dict = dict()
    
    for i in range(matrix_onlyIndex.shape[0]):
        first = matrix_onlyIndex[i]
        
        for j in range(i+1, matrix_onlyIndex.shape[0]):
            second = matrix_onlyIndex[j]
            hamming_dist = len(loci)*hamming(first, second)
            
            if hamming_dist not in hamming_dict.keys():
                hamming_dict[hamming_dist] = 1
            else:
                hamming_dict[hamming_dist] = hamming_dict[hamming_dist] + 1
                
            hammingMatrix[i,j] = 1.0*hamming_dist/len(loci)
            hammingMatrix_asdict[(ST_alleles_dict[i],ST_alleles_dict[j])] = 1.0*hamming_dist/len(loci)
            hammingMatrix[j,i] = hammingMatrix[i,j]
            hammingMatrix_asdict[(ST_alleles_dict[j],ST_alleles_dict[i])] = hammingMatrix_asdict[(ST_alleles_dict[i],ST_alleles_dict[j])]
    
    print hamming_dict
    sorted_key=sorted(hamming_dict.keys())
    plt.bar(range(len(sorted_key)), [hamming_dict[i] for i in sorted_key], align="center")
    plt.xticks(range(len(sorted_key)), sorted_key)
    plt.xlabel('Hamming Distance')
    plt.ylabel('Frequency')
    plt.title('Distribution of pairwise hamming distances \namong {} strains'.format(name))
    #plt.tight_layout()
    plt.savefig('hamDist_distribution_{}.png'.format(name), dpi=800)
    
    return hammingMatrix, hammingMatrix_asdict, matrix_onlyIndex

def split_into_clusters(link_mat,thresh,n):
   c_ts=n
   clusters={}
   for row in link_mat:
      if row[2] < thresh:
          n_1=int(row[0])
          n_2=int(row[1])

          if n_1 >= n:
             link_1=clusters[n_1]
             del(clusters[n_1])
          else:
             link_1= [n_1]

          if n_2 >= n:
             link_2=clusters[n_2]
             del(clusters[n_2])
          else:
             link_2= [n_2]

	  link_1.extend(link_2)
          clusters[c_ts] = link_1
          c_ts+=1
      else:
          return clusters

strainsAndPropFolder = "/home/glgan/Documents/Borrelia/RelevantWork/pipeline/strainsAndPropNew_23hr"
allCsv = [i for i in os.listdir(strainsAndPropFolder) if (i.startswith("SRR") and i.endswith(".csv"))]
df = pd.DataFrame()
loci = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']

for csv in allCsv:
    cols = pd.read_csv(strainsAndPropFolder+"/"+csv, nrows=1).columns
    tempDf = pd.read_csv(strainsAndPropFolder+"/"+csv, usecols=cols[1:], index_col=False)
    tempDf["Sample"] = [csv.split("_")[0]]*tempDf.shape[0]
    df = df.append(tempDf)

#df = df[df["Proportion"] > 0]
#df = df.reset_index(drop=True)

'''Number of novel and existing strains'''
df_new = df[df["New/Existing"] =="New"]
df_exist = df[df["New/Existing"] == "Existing"]
print df_exist
#print("Number of novel strains: {}".format(df_new.drop_duplicates(subset=loci).shape[0]))
#print("Number of existing strains: {}".format(df_exist.drop_duplicates(subset=loci).shape[0]))
df_justStrains=df.drop_duplicates(subset=loci)
#print("Number of strains found: {}".format(df_justStrains.shape[0]))
#df_justStrains[["ST", "New/Existing"]+ loci].to_csv("foundStrains.csv")

for ind in df_justStrains["ST"].unique():
    temp_df = df[df["ST"] == ind]

    if temp_df.shape[0] > 1:
        print temp_df

''' Generate graph where y=#of strains x=sample '''
#samp_df_new_pivot = df.replace({"New":True, "Existing":False}).pivot_table(index="Sample", columns="ST", values="New/Existing")
#samp_df_new_pivot["New strains"] = samp_df_new_pivot[samp_df_new_pivot == True].count(axis=1)
#samp_df_new_pivot["Existing strains"] = samp_df_new_pivot[samp_df_new_pivot == False].count(axis=1)
#
#samp_df_new_pivot[["Existing strains", "New strains"]].plot(kind='bar', legend=True, cmap=plt.cm.Paired)
#plt.xlabel("Sample")
#plt.ylabel("Number of strains")
#plt.title("Number of strains for each sample")
#plt.tight_layout()
#plt.savefig("sample_strainFreq.png",dpi=1000)


''' Generate graph of strain composition '''
#newStrComp = df_new.pivot(index='ST', columns='Sample', values='Proportion')
#newStrComp["Sum"] = newStrComp.sum(axis=1)
#newStrComp.sort_values(by='Sum', inplace=True, ascending=False)
#newStrComp.drop(['Sum'], axis=1).plot(kind='bar', stacked=True, legend=False, sort_columns=True, cmap=plt.cm.Paired)
#plt.xlabel('Strain type')
#plt.ylabel('Proportions')
#plt.title('Cumulated proportions of new strains in different samples')
#plt.tight_layout()
#plt.savefig("newStrainsComposition.png", dpi=1000)
#
#existStrComp = df_exist.pivot(index='ST', columns='Sample', values='Proportion')
#existStrComp["Sum"] = existStrComp.sum(axis=1)
#existStrComp.sort_values(by='Sum', inplace=True, ascending=False)
#existStrComp.drop(['Sum'], axis=1).plot(kind='bar', stacked=True, legend=False, sort_columns=True, cmap=plt.cm.Paired)
#plt.xlabel('Strain type')
#plt.ylabel('Proportions')
#plt.title('Cumulated proportions of existing strains in different samples')
#plt.tight_layout()
#plt.savefig("existingStrainsComposition.png", dpi=1000)

''' Generate new reference with half of new strains added to library  '''
#refStrains="strain_ref.txt"
#df_justStrains_new = df_new.drop_duplicates(subset=loci)[:]
#df_justStrains_new.reset_index(drop=True, inplace=True)
#writePath="new_robust"
#generateNewRefNewStr(refStrains, loci, writePath, df_justStrains_new)

''' Generate new reference with each experiment one existing strain is removed  '''
#Remove existing strain experiment
#df_justStrains_exist = df_exist.drop_duplicates(subset=loci)[:]
#df_justStrains_exist.reset_index(drop=True, inplace=True)
#for i in df_justStrains_exist.index:
#    generateNewRefRemoveExist(refStrains, df_justStrains_exist.loc[i,loci].values, loci, i, writePath)

''' Time taken for each sample for ADP '''
#timeTakenAlleleDiv = pd.read_csv("after_revision_figure/time_ilp1.csv")
#timeTakenAlleleDiv.sort_values(by="Time(sec)", inplace=True, ascending=False)
#timeTakenAlleleDiv["log(Time in sec)"] = np.log(timeTakenAlleleDiv["Time(sec)"])
#timeTakenAlleleDiv.plot(kind="bar",x="Sample", y="log(Time in sec)")
#plt.xlabel("Samples")
#plt.ylabel("Natural log of time taken in seconds")
#plt.title("Natural log of time taken for each sample when \nsolving the allele diversity problem")
#plt.tight_layout()
#plt.savefig("after_revision_figure/time_alleleDiv.png", dpi=1000)


''' Generate graph where x=ST and Y=#strains'''
#df_new = df[df["New/Existing"] =="New"]
#print df_new.drop_duplicates(subset=loci).shape
##df_new["Exist"] = [1]*df_new.shape[0]
#df_exist = df[df["New/Existing"] == "Existing"]
#df_new_pivot = df_new.pivot_table(index="ST", columns="Sample", values="Proportion")
#df_new_pivot[">1% proportion count"] = df_new_pivot[df_new_pivot > 0.01].count(axis=1)
#df_new_pivot["Count"] = df_new_pivot.count(axis=1)
#df_new_pivot["Count"] = df_new_pivot["Count"] - 1

#df_new_pivot[["Count", ">1% proportion count"]].plot(kind='bar', legend=True, cmap=plt.cm.Paired)
#plt.xlabel("Strain Type")
#plt.ylabel("Frequency")
#plt.title("Number of samples sharing each new strain")
#plt.tight_layout()
#plt.savefig("shareCount_new_beforeClust.png",dpi=1000)

#df_exist_pivot = df_exist.pivot_table(index="ST", columns="Sample", values="Proportion")
#df_exist_pivot[">1% proportion count"] = df_exist_pivot[df_exist_pivot > 0.01].count(axis=1)
#df_exist_pivot["Count"] = df_exist_pivot.count(axis=1)
#df_exist_pivot["Count"] = df_exist_pivot["Count"] - 1

#df_exist_pivot[["Count", ">1% proportion count"]].plot(kind='bar', legend=True, cmap=plt.cm.Paired)
#plt.xlabel("Strain Type")
#plt.ylabel("Frequency")
#plt.title("Number of samples sharing each existing strain")
#plt.tight_layout()
#plt.savefig("shareCount_exist.png",dpi=1000)

#When group the cherries
#df_new["Group"] = np.nan
#count=1
#df_new.loc[df_new["ST"].isin([6345, 10729]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([8344, 11611]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([16359,17609]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([14505, 16483]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([5855,6124]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([4384,12429]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([5,17896]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([5099, 5114]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([22, 224]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([261,4441,25,26,16,17652,4426,4812]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([4361,5559]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([6155,6156]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([1582,1875]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([13131,4733,12609]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([12316,18110]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([12728, 4494, 12646]), "Group"]= "G{}".format(count)
#count += 1
#df_new.loc[df_new["ST"].isin([4354, 5857]), "Group"]= "G{}".format(count)
#count += 1
#num_nan = pd.isnull(df_new["Group"]).sum()
#df_new.loc[pd.isnull(df_new["Group"]), "Group" ]= ["G{}".format(i) for i in range(count, count+num_nan)]

#grouped_df_new = df_new[["Sample", "Group"]].drop_duplicates().groupby(["Group"]).count()
#grouped_df_new.rename(columns={"Sample":"Count"}, inplace=True)
#grouped_df_new.sort_values("Count", ascending=False)["Count"].plot(kind='bar', cmap=plt.cm.Paired)
#plt.xlabel("Cluster Type")
#plt.ylabel("Frequency")
#plt.title("Number of samples sharing each cluster")
#plt.tight_layout()
#plt.savefig("shareCount_cluster.png",dpi=1000)

''' Assign group to existing strains too '''
#df_group_all = df.merge(df_new, how="outer")
#df_justStrains_exist = df_exist.drop_duplicates(subset=["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"])
#st_group = {i:"E{}".format(j) for (i,j) in itertools.izip(df_justStrains_exist["ST"].tolist(), range(1, 1+df_justStrains_exist.shape[0]))}
#df_group_all["Group"] = df_group_all["Group"].fillna(df_group_all["ST"].map(st_group))
#grouped_all_count = df_group_all[["Sample","Group"]].drop_duplicates().groupby(["Group"]).count()
#grouped_all_count.rename(columns={"Sample":"Count"}, inplace=True)
#grouped_all_count["Count"].plot(kind='bar', cmap=plt.cm.Paired)
#plt.xlabel("Group type")
#plt.ylabel("Frequency")
#plt.title("Number of samples sharing each type. G(X) represents \na cluster of new strains, E(X) represents existing strain")
#plt.tight_layout()
#plt.savefig("shareCount_all.png",dpi=1000)

#sampComparison_dict = {}
#for (i,j) in itertools.combinations(df_group_all["Sample"].unique().tolist(),2):
#    sampComparison_dict[(i,j)] = [list(set.intersection(set(df_group_all[df_group_all["Sample"] == i]["Group"].tolist()),set( df_group_all[df_group_all["Sample"] == j]["Group"].tolist())))]

#sampComparison_df = pd.DataFrame.from_dict(sampComparison_dict, orient='index')
#sampComparison_df = sampComparison_df.rename(columns={0:"Common"})
#sampComparison_df["Count"] = sampComparison_df.apply(lambda x: len(x["Common"]), axis=1)
#sampComparison_df_nonzero = sampComparison_df[sampComparison_df["Count"] != 0].reset_index()
''' Generate graph where y=#of strains x=sample '''
#samp_df_new_pivot = df.replace({"New":True, "Existing":False}).pivot_table(index="Sample", columns="ST", values="New/Existing")
#samp_df_new_pivot["New strains"] = samp_df_new_pivot[samp_df_new_pivot == True].count(axis=1)
#samp_df_new_pivot["Existing strains"] = samp_df_new_pivot[samp_df_new_pivot == False].count(axis=1)

#samp_df_new_pivot[["Existing strains", "New strains"]].plot(kind='bar', legend=True, cmap=plt.cm.Paired)
#plt.xlabel("Sample")
#plt.ylabel("Number of strains")
#plt.title("Number of strains for each sample")
#plt.tight_layout()
#plt.savefig("sample_strainFreq.png",dpi=1000)
#sample_group_df_new = pd.pivot_table(df_new, index="Sample", columns="Group", values="Exist")
#sample_group_df_new.fillna(0, inplace=True)
#sample_group_df_new['compressed'] = sample_group_df_new.apply(lambda x: ''.join([ str(v) for v in x ]),1)
#sample_setOfStr = sample_group_df_new.groupby('compressed').apply(lambda x: x.index.tolist())
#sample_setOfStr = sample_setOfStr.reset_index()
#sample_setOfStr[0].to_csv("groupShared.csv")
#df_new = df_new[df_new["Proportion"] > 0]
#df_new.sort(columns=["Sample"], inplace=True)
#df_exist = df_exist[df_exist["Proportion"] > 0]
#df_exist.sort(columns=["Sample"], inplace=True)
#df_new.to_csv("newSamples.csv")
#Just consider strain types
#df_justStrains_new = df_new.drop_duplicates(subset=["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"])[:]
#Any samples sharing strains
#df_diff_new = df_new[~df_new.index.isin(df_justStrains_new.index)]
#df_justStrains_exist = df_exist.drop_duplicates(subset=["clpA", "clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA"])
#df_diff_exist = df_exist[~df_exist.index.isin(df_justStrains_exist.index)]
#df_justStrains[["New/Existing", "clpA","clpX", "nifS", "pepX", "pyrG", "recG", "rplB", "uvrA" ]].to_csv("newStrains.csv")

##Remove existing strain experiment
#for i in df_justStrains_exist.index:
#    generateNewRefRemoveExist(refStrains, df_justStrains_exist.loc[i,loci].values, loci, i, "/home/glgan/Documents/Borrelia/pipeline")

#newStrComp = df_new.pivot(index='ST', columns='Sample', values='Proportion')
#newStrComp["Sum"] = newStrComp.sum(axis=1)
#newStrComp.sort_values(by='Sum', inplace=True, ascending=False)
#newStrComp.drop(['Sum'], axis=1).plot(kind='bar', stacked=True, legend=False, sort_columns=True, cmap=plt.cm.Paired)
#plt.xlabel('Strain type')
#plt.ylabel('Proportions')
#plt.title('Cumulated proportions of new strains in different samples')
#plt.tight_layout()
#plt.savefig("newStrainsComposition_err25.png", dpi=1000)

#existStrComp = df_exist.pivot(index='ST', columns='Sample', values='Proportion')
#existStrComp["Sum"] = existStrComp.sum(axis=1)
#existStrComp.sort_values(by='Sum', inplace=True, ascending=False)
#existStrComp.drop(['Sum'], axis=1).plot(kind='bar', stacked=True, legend=False, sort_columns=True, cmap=plt.cm.Paired)
#plt.xlabel('Strain type')
#plt.ylabel('Proportions')
#plt.title('Cumulated proportions of existing strains in different samples')
#plt.tight_layout()
#plt.savefig("existingStrainsComposition_err25.png", dpi=1000)

