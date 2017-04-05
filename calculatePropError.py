#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 18:32:04 2017

@author: stanleygan
"""
import os
import glob
import re
import numpy as np
import pandas as pd

def calculateRMSError(predictedFile, trueFile):
    pred = pd.read_csv(predictedFile, names=["Variant", "Predicted Proportion"])
    true = pd.read_csv(trueFile, names=["Variant", "True Proportion"])
    
    compare = pred.merge(true, indicator=True, how="outer")
    
    common = compare["_merge"] == "both"
    commonDF = compare[common]
    diff = np.asarray( (commonDF["Predicted Proportion"] - commonDF["True Proportion"]).tolist() )
    rmse = np.sqrt( np.mean(np.square(diff)) )
    print("For those variants which are predicted correctly, the Root Mean Square Error of the proportions is:{}".format(rmse))
    
    predictedButFalse = compare["_merge"] == "left_only"
    trueButNotPred = compare["_merge"] == "right_only"
    
    if not all(~predictedButFalse):
        print("Following variants are predicted but are not true variants:")
        print(compare[predictedButFalse]["Variant"].tolist())
    
    if not all(~trueButNotPred):
        print("Following variants are true but not predicted:")
        print(compare[trueButNotPred]["Variant"].tolist())
    
    print("\n")

predictedDir = "/home/stanleygan/Documents/Borrelia/data/predictedVar/"
trueDir = "/home/stanleygan/Documents/Borrelia/data/simData/trueVar/"

os.chdir(predictedDir)
predictedCsv = [f for f in glob.glob("*proportions.csv")]
os.chdir(trueDir)
trueCsv = [f for f in glob.glob("*proportions.csv")]

temp = re.compile("[a-zA-Z]{4}_[0-9]*")
match = [temp.findall(s) for s in trueCsv]

for name in match:
    predictedFile = [re.compile(name[0]).findall(s) for s in predictedCsv if len(re.compile(name[0]).findall(s)) != 0]
    trueFile =  [re.compile(name[0]).findall(s) for s in trueCsv if len(re.compile(name[0]).findall(s)) != 0]
    
    print("For " + name[0] + ":")
    calculateRMSError(predictedDir+predictedFile[0][0]+"_predicted_proportions.csv", trueDir+trueFile[0][0]+"_proportions.csv")
    