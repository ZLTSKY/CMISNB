import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.stats as stat
from pyitlib import discrete_random_variable as drv
import multiprocessing
import math
import knncmi as k

    
begin=time.asctime()
print("Begin time: "+begin)

param={}
for i in range(1,len(sys.argv)):
    t=sys.argv[i].split("=")
    param[t[0].lower()]=t[1]

help_msg="""
usage: python knncmi_network.py -process=process_value -background=background_network_file -ref=reference_sample_file  -sample=sample_data_file -out=results_output_fold
Options and arguments:
-proportion: the genes whose expression value was 0 in more than proportion samples/cells
-process: namesumber of work processes used
-background : background network to calculate the deltaCMI of edges based on the network
-ref : the expression profile of reference samples
-sample : the expression profile for the sample to be constructed the CMISNB
-out : the directory to store the CMISNB
"""
if "-help" in param.keys() or "-h" in param.keys():
    print(help_msg)


if "-ref" not in param.keys() or "-sample" not in param.keys() or "-background" not in param.keys() or "-process" not in param.keys():
    print("Parameter missing!")
    print(help_msg)
    exit()
reference_file=param["-ref"]
sample_file=param["-sample"]
background=param["-background"]
process=int(param["-process"])


if "-out" not in param.keys():    
    fold="."
else:
    fold=param["-out"]

if not os.path.exists(fold):
    os.mkdir(fold)


###############################################################################################################
#############################                                          ########################################
#############################                  cmi代码                 ########################################
#############################                                          ########################################
#############################        CMI代码所需的库主要为knncmi       ########################################
#############################                                          ########################################
###############################################################################################################


def cmi(X,Y,Z):
    if len(Z)>0:
        W=[X]+[Y]+[Z]
        W = np.transpose(W)
        A = pd.DataFrame(W)
        CMI = k.cmi([0], [1], [2],3, A)
    else:
        W=[X]+[Y]
        W = np.transpose(W)
        A = pd.DataFrame(W)
        CMI = k.cmi([0], [1], [],3, A)      
    return CMI

###############################################################################################################
#############################                                          ########################################
#############################                  cmi代码                 ########################################
#############################                                          ########################################
###############################################################################################################

def parallel_procedure(index,names,gene_filter,gene_correspondence,data,ref,ref_count):
    print("inner process",index)
    fw=open(fold +os.sep+"cmi_"+ names[index]+".txt","w")
    fw.write("Gene1\tGene2\tmean_deltaCMI\tCMI\n")
    fw.close()
    for i in range(len(gene_filter)):
        s=gene_filter[i].split("_")
        deltacmi=[]
        if len(gene_correspondence[gene_filter[i]])>0:
            for j in gene_correspondence[gene_filter[i]]:
                rc = cmi(ref[s[0]],ref[s[1]],ref[j])
                rc1 = cmi(ref[s[0]]+[data[s[0]][index]],ref[s[1]]+[data[s[1]][index]],ref[j]+[data[j][index]])-rc
                deltacmi.append(rc1)
        else:
            rc = cmi(ref[s[0]],ref[s[1]],[])
            rc1 = cmi(ref[s[0]]+[data[s[0]][index]],ref[s[1]]+[data[s[1]][index]],[])-rc
            deltacmi.append(rc1)
        mean_deltacmi=np.mean(deltacmi)
        fw=open(fold +os.sep+"cmi_"+ names[index]+".txt","a")
        fw.write(s[0]+"\t"+s[1]+"\t"+str(mean_deltacmi))
        for num in range(len(deltacmi)):
            fw.write("\t"+str(deltacmi[num]))
        fw.write('\n')
        fw.close()



if __name__=="__main__":
    ref={}
    ref_count=0
    f=open(reference_file)
    flag=0
    for p in f:
        flag+=1
        t=p.split()
        if flag==1:
            ref_count=len(t)
            continue
        if t[0] not in ref.keys():
            ref[t[0]]=[float(t[i]) for i in range(1,len(t))]
        else:
            print("Error in ",t[0])
    f.close()
    genes=list(ref.keys())

    data={}
    f=open(sample_file)
    flag=0
    for p in f:
        flag+=1
        t=p.split()
        if flag==1:
            names=t[0:]
            continue
        if t[0] not in data.keys():
            data[t[0]]=[float(t[i]) for i in range(1,len(t))]
        else:
            print("Error in ",t[0])
    f.close()


    gene_correspondence = {}
    f=open(background)
    flag=0
    for p in f:
        t = p.split()
        if t[0] not in gene_correspondence.keys():
            gene_correspondence[t[0]]=[t[i] for i in range(1,len(t))]
        else:
            print("Error in",t[0])
    f.close()
    gene_filter=list(gene_correspondence.keys())

    pool=multiprocessing.Pool(process)
    for index in range(len(names)):#左闭右开
        pool.apply_async(parallel_procedure,(index,names,gene_filter,gene_correspondence,data,ref,ref_count,))
    print('Waiting for all subprocesses done...')  
    pool.close()
    pool.join()




print("Begin time: "+begin)
print("End time: "+time.asctime())
