import numpy as np
import pandas as pd
import os
import igraph as ig

dirPath=r'xxxxxxxxxxxxxxxx' # The address of the input file
storagePath=r'xxxxxxxxxxxxxxxx' # The address of the output file
background_file='xxxxxxxxxxxx' # After screening, the background network that meets the expression matrix and comprehensive score should be careful to remove repeated edges, such as A-B,B-A
gene_partial_background='xxxxxxxxxxxxx' # The name of output file

g= ig.Graph.Read_Ncol(dirPath + os.sep + background_file,names=True,directed=False)

f=open(dirPath + os.sep + background_file,'r')
new_ppi=[]
for p in f:
    t = p.split()
    new_ppi.append(t)
f.close()

fw=open(storagePath + os.sep + gene_partial_background,'w')
fw.close()
for number in range(len(new_ppi)):
    print(new_ppi[number])
    gene1_neighbors=g.neighbors(new_ppi[number][0])#
    gene2_neighbors=g.neighbors(new_ppi[number][1])#
    gene12_neighbors=[]
    [gene12_neighbors.append(x) for x in gene1_neighbors if x in gene2_neighbors if x not in gene12_neighbors]
    gene=[]
    for i_z in range(len(gene12_neighbors)):
        gene.append(g.vs[gene12_neighbors[i_z]]['name'])
    gene_partial=" ".join(gene)
    fw=open(storagePath + os.sep + gene_partial_background,'a')
    fw.write(str(new_ppi[number][0])+"_"+str(new_ppi[number][1])+ "\t" + "\t" +gene_partial + "\n")
    fw.close()