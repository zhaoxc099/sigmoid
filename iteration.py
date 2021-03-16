from sympy import symarray
import numpy as np
import xlrd
import math
import time
import pandas as pd
import csv

start_time = time.time()

with open('data.csv','r') as csvfile:
    reader = csv.reader(csvfile)
    rows= [col for col in reader]

rows=np.array(rows)
sigma=rows.astype(np.float)


print("--- %s seconds ---" % (time.time() - start_time))

M=len(sigma)         #number of rows
N=len(sigma[0])    #number of columns

print(M)          #M=59944
print(N)           #N=158 

iter           = 400000
eta            = 0.0001

K=30        #need to change from 5~50
#M=50000
sigma=sigma[0:M]
Beta           = 0.1*np.random.rand(M,K)
Energy         = 0.1*np.random.rand(K,N)


for z in range(iter):
    
    theta  = np.dot(Beta,Energy)
    pi = 1/(1+np.exp(theta))
    delta = pi-sigma
    


    init_E=np.dot(Beta.T,delta)
    init_B=np.dot(Energy,delta.T)
    
    
    Beta=Beta+(eta*init_B).T
    Energy=Energy+eta*init_E

    if z % 1000 == 0:
        p=0
        for i in range(M):
            for j in range(N):
                if pi[i][j]>0.999999999:
                    p=p+sigma[i][j]*math.log(pi[i][j])
                elif pi[i][j]<0.000000001:
                    p=p+(1-sigma[i][j])*math.log(1-pi[i][j])
                else:
                    p=p+sigma[i][j]*math.log(pi[i][j])+(1-sigma[i][j])*math.log(1-pi[i][j])
        print("iteration=%s"%z)
        print("converge=%s"%p)  
        print("--- %s seconds ---" % (time.time() - start_time))
        
        e=0
        for i in Energy:
            for j in i:
                e=e+j*j
        de=0
        for i in init_E:
            for j in i:
                de=de+j*j
        print((de/e)**0.5)
        
        if (de/e)**0.5<0.001:
            break
    
        np.savetxt("dBk5.csv", init_B, delimiter=',')
        np.savetxt("dEk5.csv", init_E, delimiter=',')
        np.savetxt("Bk5.csv", Beta, delimiter=',')
        np.savetxt("Ek5.csv", Energy, delimiter=',')
        np.savetxt("pik5.csv", pi, delimiter=',')






