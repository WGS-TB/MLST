import pandas as pd
import numpy as np
import os
import networkx as nx
import matplotlib.pyplot as plt

def returnNameToStrainDict(df, loci):
    assert('newNaming' in df.columns)
    
    d = dict()
    for row in df[loci+['newNaming']].itertuples(index=False):
        d[list(row)[-1]] = list(row)[:-1]

    return d       
 
def returnDistMatrix(d,loci):
    distMat = np.zeros((len(d),len(d)))

    for row, st1 in enumerate(d.keys()):
        for col, st2 in enumerate(d.keys()):
            if st1 == st2:
                distMat[row, col] = 0
            else:
                distMat[row,col] = computeStrainDist(d[st1], d[st2],loci)

    return distMat

def computeStrainDist(st1, st2,loci):
    dist = 0
    for i,gene in enumerate(loci):
        dist += computeGeneDist(gene,st1[i], st2[i])  

    return dist

def computeGeneDist(gene, a1, a2):
    distanceMat = genes_distance[gene]         
    return int(distanceMat.loc[distanceMat['level_0'] == a1][a2])
 
project_path = os.getcwd()
strainsAndPropFolder = "/home/glgan/Documents/Borrelia/RelevantWork/pipeline/newS"
allCsv = [i for i in os.listdir(strainsAndPropFolder) if (i.startswith("SRR") and i.endswith(".csv"))]
df = pd.DataFrame()
loci = ['clpA', 'clpX', 'nifS', 'pepX', 'pyrG', 'recG', 'rplB', 'uvrA']


#editDistanceMatrices for each gene
clpA_df = pd.read_csv(os.path.join(project_path,'editDist','editDistanceMatrix_clpA.csv'), sep=",")
#clpX
clpX_df = pd.read_csv(os.path.join(project_path, 'editDist','editDistanceMatrix_clpX.csv'), sep=",")
#nifS
nifS_df = pd.read_csv(os.path.join(project_path, 'editDist','editDistanceMatrix_nifS.csv'), sep=",")
#pepX
pepX_df = pd.read_csv(os.path.join(project_path, 'editDist','editDistanceMatrix_pepX.csv'), sep=",")
#pyrG
pyrG_df = pd.read_csv(os.path.join(project_path, 'editDist','editDistanceMatrix_pyrG.csv'), sep=",")
#recG
recG_df = pd.read_csv(os.path.join(project_path, 'editDist','editDistanceMatrix_recG.csv'), sep=",")
#rplB
rplB_df = pd.read_csv(os.path.join(project_path, 'editDist','editDistanceMatrix_rplB.csv'), sep=",")
#uvrA
uvrA_df = pd.read_csv(os.path.join(project_path, 'editDist','editDistanceMatrix_uvrA.csv'), sep=",")
#dictionary
genes_distance = {'clpA':clpA_df, 'clpX':clpX_df, 'nifS':nifS_df, 'pepX':pepX_df, 'pyrG':pyrG_df, 'recG':recG_df, 'rplB':rplB_df, 'uvrA':uvrA_df}

for csv in allCsv:
    cols = pd.read_csv(strainsAndPropFolder+"/"+csv, nrows=1).columns
    tempDf = pd.read_csv(strainsAndPropFolder+"/"+csv, usecols=cols[1:], index_col=False)
    tempDf["Sample"] = [csv.split("_")[0]]*tempDf.shape[0]
    df = df.append(tempDf)

df = df.drop_duplicates(subset=loci).reset_index()

newNaming = list()
nIdx = 0
eIdx = 0

for nOrE in df['New/Existing'].tolist():
    if nOrE == 'New':
        newNaming.append( "n{}".format(nIdx))
        nIdx += 1
    elif nOrE == 'Existing':
        newNaming.append( "e{}".format(eIdx))
        eIdx += 1

df['newNaming'] = newNaming

nameDict = returnNameToStrainDict(df, loci)    
distMat = returnDistMatrix(nameDict, loci)
distMat_df = pd.DataFrame(distMat, index=nameDict.keys(), columns=nameDict.keys())
graph = nx.from_pandas_adjacency(distMat_df)
mst = nx.minimum_spanning_tree(graph)
#node_pos = nx.get_node_attributes(mst,'pos')
node_pos = nx.spring_layout(mst, k=4*1/np.sqrt(len(mst.nodes())), iterations=80, random_state=1992)
labels = nx.get_edge_attributes(mst, 'weight')
plt.figure()
nx.draw_networkx_edge_labels(mst, node_pos, edge_labels=labels, font_size=4)
nx.draw_networkx(mst,pos=node_pos,node_size=10,font_size=1, alpha=0.65)
#nx.draw_networkx(mst,pos=node_pos,node_size=90,font_size=5, alpha=0.65)
plt.savefig("MST.png", dpi=350)
