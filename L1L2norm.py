import numpy as np
import math
import time
import csv


def l2norm(K):
    with open('Bk%s.csv' %K,'r') as csvfile:
        reader = csv.reader(csvfile)
        rows= [col for col in reader]
    
    rows=np.array(rows)
    Beta=rows.astype(np.float)
    
    with open('Ek%s.csv'%K,'r') as csvfile:
        reader = csv.reader(csvfile)
        rows= [col for col in reader]
    
    rows=np.array(rows)
    Energy=rows.astype(np.float)
    
    Beta1=[]
    Energy1=[]
    
    for i in range(K):
        n=np.linalg.norm(Beta.T[i])
        Beta1.append(Beta.T[i]/n)
        Energy1.append(Energy[i]*n)
        
    Beta=np.array(Beta1).T
    Energy=np.array(Energy1)
        
    np.savetxt("L2Bk%s.csv" %K, Beta, delimiter=',')
    np.savetxt("L2Ek%s.csv" %K, Energy, delimiter=',')
    return


def l1norm(K):
    with open('Bk%s.csv' %K,'r') as csvfile:
        reader = csv.reader(csvfile)
        rows= [col for col in reader]
    
    rows=np.array(rows)
    Beta=rows.astype(np.float)
    
    with open('Ek%s.csv'%K,'r') as csvfile:
        reader = csv.reader(csvfile)
        rows= [col for col in reader]
    
    rows=np.array(rows)
    Energy=rows.astype(np.float)
    
    Beta1=[]
    Energy1=[]
    
    for i in range(K):
        n=np.linalg.norm(Beta.T[i],ord=1)
        Beta1.append(Beta.T[i]/n)
        Energy1.append(Energy[i]*n)
        
    Beta=np.array(Beta1).T
    Energy=np.array(Energy1)
        
    np.savetxt("L1Bk%s.csv" %K, Beta, delimiter=',')
    np.savetxt("L1Ek%s.csv" %K, Energy, delimiter=',')
    return


for i in [5,10,15,20,25,30,35,40,45]:
    l2norm(i)
    l1norm(i)
    










