#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 17:02:13 2017

@author: glgan
"""

import pandas as pd
import os
import matplotlib.pyplot as plt

def returnStrainsAndPropDF(folder):
    strainsAndPropFolder = folder
    allCsv = [i for i in os.listdir(strainsAndPropFolder) if (i.startswith("SRR") and i.endswith(".csv"))]
    df = pd.DataFrame()
    
    for csv in allCsv:
        cols = pd.read_csv(strainsAndPropFolder+"/"+csv, nrows=1).columns
        tempDf = pd.read_csv(strainsAndPropFolder+"/"+csv, usecols=cols[1:], index_col=False)
        tempDf["Sample"] = [csv.split("_")[0]]*tempDf.shape[0]
        df = df.append(tempDf)
        
    df = df[df["Proportion"] > 0]
    df = df.reset_index(drop=True)
    
    return df

def booleanAllAddedFound(added, current):
    df = returnStrainsAndPropDF(current)
    df_added = pd.read_csv(added,usecols=range(1,len(loci)+1))
    for name in loci:
        df_added["%s" %name] = name + "_" + df_added["%s" %name].astype(str)
    merged = df_added[loci].merge(df.drop_duplicates(subset=loci), indicator=True, how='left')
    if df_added.shape[0] == (merged["_merge"] == "both").sum():
        return True
    else:
        return False

def prec_recall(df_ori, df_ran):
    df_merged = df_ori.merge(df_ran, on=loci, indicator=True, how='outer')[loci+["_merge", "New/Existing_x", "New/Existing_y"]].drop_duplicates(subset=loci)
    tp = (df_merged["_merge"] == "both").sum()
    precision = 1.0*tp/(tp+ (df_merged["_merge"] == "right_only").sum())
    recall = 1.0*tp/(tp+ (df_merged["_merge"] == "left_only").sum())
    
    return precision, recall

def singleObs_tvd(obs):
    if obs["_merge"] == "both":
        return abs(obs["Proportion_x"] - obs["Proportion_y"])
    elif obs["_merge"] == "left_only":
        return obs["Proportion_x"]
    else:
        return obs["Proportion_y"]

def totalVarDist(df_ori, df_ran):
    df_merged = df_ori.merge(df_ran, on=loci+["Sample"], indicator=True, how='outer')
    df_merged["Diff"] = df_merged.apply(singleObs_tvd, axis=1)
    tvd = df_merged["Diff"].sum()/2.0
    
    return tvd
    
loci = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']
currentPath = "/home/glgan/Documents/Borrelia/pipeline"
''' ===== Analysis ===='''
df_ori = returnStrainsAndPropDF("/home/glgan/Documents/Borrelia/pipeline/strainsAndProp_3hr")

''' All removed new strains found? '''
for i in range(1,11):
    print(booleanAllAddedFound(currentPath+"/added_to_ref{}.csv".format(i), currentPath+"/robust_new/strainsAndProp_3hr_ran{}".format(i)))

''' Removed existing strain found? '''
df_remExist = returnStrainsAndPropDF(currentPath+"/robust_exist/e_strainsAndProp_3hr_ref{}".format(80))
merged_exist = df_justStrains_exist[loci].merge(df_remExist.drop_duplicates(subset=loci), indicator=True, how='left')
print merged_exist

''' Precision, recall and tvd? for first experiment'''
precision = list()
recall = list()
tvd = list()

for i in range(1,11):
    df_exp = returnStrainsAndPropDF(currentPath+"/robust_new/strainsAndProp_3hr_ran{}".format(i))
    p, r = prec_recall(df_ori, df_exp)
    t = totalVarDist(df_ori, df_exp)
    precision.append(p)
    recall.append(r)
    tvd.append(t)
    
plt.figure()
plt.boxplot([precision, recall])
plt.xticks([1,2], ["Precision", "Recall"])
plt.xlabel("Statistics")
plt.ylabel("Value")
plt.title("Precision and recall across 10 experiments \nwhere half of new strains are added")
plt.tight_layout()
plt.savefig("prec_recall_addNewStr.png", dpi=1000)

plt.figure()
plt.boxplot(tvd)
plt.xticks([1], ["Total Variation Distance"])
plt.ylabel("Value")
plt.title("Total Variation Distance across 10 experiments \nwhere half of new strains are added")
plt.tight_layout()
plt.savefig("tvd_addNewStr.png", dpi=1000)
    
''' Precision...for 2nd experiment '''
indexRem = [0,1,3,4,11,23,26,42,45,47,49,50,62,80]
precision_2 = list()
recall_2 = list()
tvd_2 = list()
for i in indexRem:
    df_remExist = returnStrainsAndPropDF(currentPath+"/robust_exist/e_strainsAndProp_3hr_ref{}".format(i))
    p,r = prec_recall(df_ori, df_remExist)
    t = totalVarDist(df_ori, df_remExist)
    precision_2.append(p)
    recall_2.append(r)
    tvd_2.append(t)
    
plt.figure()
plt.boxplot([precision_2, recall_2])
plt.xticks([1,2], ["Precision", "Recall"])
plt.xlabel("Statistics")
plt.ylabel("Value")
plt.title("Precision and recall across 14 experiments \nwhere an existing strain is removed")
plt.tight_layout()
plt.savefig("prec_recall_removeExisting.png", dpi=1000)

plt.figure()
plt.boxplot(tvd_2)
plt.xticks([1], ["Total Variation Distance"])
plt.ylabel("Value")
plt.title("Total Variation Distance across 10 experiments \nwhere an existing strain is removed")
plt.tight_layout()
plt.savefig("tvd_removeExisting.png", dpi=1000)

''' Compare samples before and after '''
#i=1
#df_ran1 = returnStrainsAndPropDF(currentPath+"/robust_new/strainsAndProp_3hr_ran{}".format(i))
#df_merged = df_ori.merge(df_ran1, on=loci, indicator=True, how='outer')
