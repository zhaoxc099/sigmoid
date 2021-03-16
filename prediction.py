import numpy as np
import xlrd
import math
import time
import csv
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch
import scipy
import scipy.io
from sklearn.metrics import jaccard_score

start_time = time.time()

# load csv file----------------------------------------------------------------

def openfile(filename):
    with open(filename,'r') as csvfile:
        reader = csv.reader(csvfile)
        rows= [col for col in reader]
    rows=np.array(rows)
    rows=rows.astype(np.float)  #raw data
    return rows

#-----------------------------------------------------------------------
K=30 # k value
X=openfile('L2Ek%s.csv'%K) #load normlized E matrix


#----------------------------------------

linked = linkage(X.T, 'ward')         # hierarchical clustering

#---------------------------------------------------------------------------

#tune clustering with parameter y

def find_idx(y):
    idx=sch.fcluster(linked, y*X.max(),'distance').tolist()
    R=max(idx)
    
    index=scipy.zeros([R,0]).tolist()
    for n in range(N):
        tmp=idx[n]
        index[tmp-1].append(n)
    return index


# load testing set---------------------------------------------------------------
sigma=openfile("data_test.csv")
M,N=sigma.shape



#split testing set into known (rxn_80) and unkown groups (rxn_20)-------------

def split(m,p,sigma,N):
    rxn=[]
    non_rxn=[]
    for n in range(N):
        if sigma[m][n]==1:
            rxn.append(n) 
        else:
            non_rxn.append(n)
            
    rxn_80=[]
    rxn_20=[]
    for r in rxn:
        if random.random()<p:
            rxn_80.append(r)
        else:
            rxn_20.append(r)
    return rxn_80, rxn_20


# predict unkown reactions based on known reactions, and evaluate prediction with TPR and FPR


def roc_hierachical(m,rxn_80,rxn_20,sigma,index,N):
    
    FP=[]
    TP=[]
    FN=[]
    TN=[]

    predict=list(np.zeros(N))
    L=len(rxn_80)
    Q=len(rxn_20)

    for i in index:
        if len(set(i)&set(rxn_80))>=1:
            for k in i:
                predict[k]=1


    for n in range(N):
        if n in rxn_80:
            pass
        elif (predict[n]==1 and sigma[m][n]==1):
            TP.append(n)
        elif (predict[n]==1 and sigma[m][n]==0):
            FP.append(n)
        elif (predict[n]==0 and sigma[m][n]==1):
            FN.append(n)
        elif (predict[n]==0 and sigma[m][n]==0):
            TN.append(n)           

#    PRC=len(TP)/(len(TP)+len(FP))
    TPR=len(TP)/(len(TP)+len(FN))
    FPR=len(FP)/(len(TN)+len(FP))


    return TPR, FPR, predict



# An example to set up parameters

p=0.9 # 90% reactions are known
y=0.1 # change this number to tune clustering level
index=find_idx(y)
TPR=[]
FPR=[]
Predict_matrix=[]

for m in range(M):                #loop all testing samples
    rxn_80, rxn_20= split(m,p,sigma,N)
    tmp=roc_hierachical(m,rxn_80,rxn_20,sigma,index,N)
    TPR.append(tmp[0])
    FPR.append(tmp[1])
    Predict_matrix.append(tmp[2])

