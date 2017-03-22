# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 17:16:05 2017

@author: Nikhil Krishnan
"""

# Modified by Nara Yoon on Mar. 16, 2017

import numpy as np
#from copy import copy
#from copy import deepcopy
#from random import random
#from random import choice
#from random import randint
#import matplotlib.pyplot as plt
import pandas

from scipy.stats import rankdata    


AMP = [1.851, 2.082, 1.948, 2.434, 2.024, 2.198, 2.033, 0.034, 1.57, 2.165, 0.051, 0.083, 2.186, 2.322, 0.088, 2.821]
AM = [1.778, 1.782, 2.042, 1.752, 1.448, 1.544, 1.184, 0.063, 1.72, 2.008, 1.799, 2.005, 1.557, 2.247, 1.768, 2.047]
CEC = [2.258, 1.996, 2.151, 2.648, 2.396, 1.846, 2.23, 0.214, 0.234, 0.172, 2.242, 0.093, 2.15, 0.095, 2.64, 0.516]
CTX = [0.16, 0.085, 1.936, 2.348, 1.653, 0.138, 2.295, 2.269, 0.185, 0.14, 1.969, 0.203, 0.225, 0.092, 0.119, 2.412]
ZOX = [0.993, 0.805, 2.069, 2.683, 1.698, 2.01, 2.138, 2.688, 1.106, 1.171, 1.894, 0.681, 1.116, 1.105, 1.103, 2.591]
CXM = [1.748, 1.7, 2.07, 1.938, 2.94, 2.173, 2.918, 3.272, 0.423, 1.578, 1.911, 2.754, 2.024, 1.678, 1.591, 2.923]
CRO = [1.092, 0.287, 2.554, 3.042, 2.88, 0.656, 2.732, 0.436, 0.83, 0.54, 3.173, 1.153, 1.407, 0.751, 2.74, 3.227]
AMC = [1.435, 1.573, 1.061, 1.457, 1.672, 1.625, 0.073, 0.068, 1.417, 1.351, 1.538, 1.59, 1.377, 1.914, 1.307, 1.728]
CAZ = [2.134, 2.656, 2.618, 2.688, 2.042, 2.756, 2.924, 0.251, 0.288, 0.576, 1.604, 1.378, 2.63, 2.677, 2.893, 2.563]
CTT = [2.125, 1.922, 2.804, 0.588, 3.291, 2.888, 3.082, 3.508, 3.238, 2.966, 2.883, 0.89, 0.546, 3.181, 3.193, 2.543]
SAM = [1.879, 2.533, 0.133, 0.094, 2.456, 2.437, 0.083, 0.094, 2.198, 2.57, 2.308, 2.886, 2.504, 3.002, 2.528, 3.453]
CPR = [1.743, 1.662, 1.763, 1.785, 2.018, 2.05, 2.042, 0.218, 1.553, 0.256, 0.165, 0.221, 0.223, 0.239, 1.811, 0.288]
CPD = [0.595, 0.245, 2.604, 3.043, 1.761, 1.471, 2.91, 3.096, 0.432, 0.388, 2.651, 1.103, 0.638, 0.986, 0.963, 3.268]
TZP = [2.679, 2.906, 2.427, 0.141, 3.038, 3.309, 2.528, 0.143, 2.709, 2.5, 0.172, 0.093, 2.453, 2.739, 0.609, 0.171]
FEP = [2.59, 2.572, 2.393, 2.832, 2.44, 2.808, 2.652, 0.611, 2.067, 2.446, 2.957, 2.633, 2.735, 2.863, 2.796, 3.203]
   
ALL = [AMP, AM, CEC, CTX, ZOX, CXM, CRO, AMC, CAZ, CTT, SAM, CPR, CPD, TZP, FEP]

#ALLreshape = np.reshape(ALLreshape,(15,2,2,2,2))
ALLreshape = np.reshape(ALL,(15,2,2,2,2))    #-------changed



#-------------------------------------------------------------------------------- changed




#-------------------------------------------------------------------------------- newly included function
# function to classify: magnitude / simple sign / reciprocal sign
def EpisClass(array22):
    ArrayFlat = array22.flatten()
    ArraySort = rankdata(ArrayFlat)
    Diag1 = np.sort(ArraySort[1:3]).tolist()
    Diag2 = np.sort(ArraySort[0:4:3]).tolist()
    
    #print(ArrayFlat,ArraySort,Diag1,Diag2)
    
    if Diag1==[3,4] or Diag2==[3,4]:
        a = 0
    elif Diag1==[1,4] or Diag2==[2,3]:
        a = 1
    else:
        a = 2
    ## 0 = reciprocal, 1 = magnitude, #2 = simple sign       
    return a

a = np.sum(np.transpose(ALLreshape,(0,1,2,3,4)).reshape((15,2,2,4)),axis=3)
b = np.sum(np.transpose(ALLreshape,(0,1,3,2,4)).reshape((15,2,2,4)),axis=3)
c = np.sum(np.transpose(ALLreshape,(0,1,4,2,3)).reshape((15,2,2,4)),axis=3)
d = np.sum(np.transpose(ALLreshape,(0,2,3,1,4)).reshape((15,2,2,4)),axis=3)
e = np.sum(np.transpose(ALLreshape,(0,2,4,1,3)).reshape((15,2,2,4)),axis=3)
f = np.sum(np.transpose(ALLreshape,(0,3,4,1,2)).reshape((15,2,2,4)),axis=3)

    
EpiSum = np.zeros((6,15)) ##results for epistatic categorization will be stored in this array
trans = ['a','b','c','d','e','f']
for j in range(0,len(trans)):    
    for i in range(0,15):
        EpiSum[j,i] = EpisClass(eval(trans[j])[i])
#print(EpiSum)

             

a = np.sum(np.exp(np.transpose(ALLreshape,(0,1,2,3,4)).reshape((15,2,2,4))),axis=3)
b = np.sum(np.exp(np.transpose(ALLreshape,(0,1,3,2,4)).reshape((15,2,2,4))),axis=3)
c = np.sum(np.exp(np.transpose(ALLreshape,(0,1,4,2,3)).reshape((15,2,2,4))),axis=3)
d = np.sum(np.exp(np.transpose(ALLreshape,(0,2,3,1,4)).reshape((15,2,2,4))),axis=3)
e = np.sum(np.exp(np.transpose(ALLreshape,(0,2,4,1,3)).reshape((15,2,2,4))),axis=3)
f = np.sum(np.exp(np.transpose(ALLreshape,(0,3,4,1,2)).reshape((15,2,2,4))),axis=3)
        
EpiExp = np.zeros((6,15)) ##results for epistatic categorization will be stored in this array
trans = ['a','b','c','d','e','f']
for j in range(0,len(trans)):    
    for i in range(0,15):
        EpiExp[j,i] = EpisClass(eval(trans[j])[i])
#print(EpiExp)



    
episclasses = ['reciprocal','magnitude','simple']
transpositions = ['1,2,3,4','1,3,2,4','1,4,2,3','2,3,1,4','2,4,1,3','3,4,1,2']
sumdata = [[EpiSum[0].tolist().count(0),EpiSum[0].tolist().count(1),EpiSum[0].tolist().count(2)], [EpiSum[1].tolist().count(0),EpiSum[1].tolist().count(1),EpiSum[1].tolist().count(2)],[EpiSum[2].tolist().count(0),EpiSum[2].tolist().count(1),EpiSum[2].tolist().count(2)],[EpiSum[3].tolist().count(0),EpiSum[3].tolist().count(1),EpiSum[3].tolist().count(2)],[EpiSum[4].tolist().count(0),EpiSum[4].tolist().count(1),EpiSum[4].tolist().count(2)],[EpiSum[5].tolist().count(0),EpiSum[5].tolist().count(1),EpiSum[5].tolist().count(2)]]
expdata = [[EpiExp[0].tolist().count(0),EpiExp[0].tolist().count(1),EpiExp[0].tolist().count(2)], [EpiExp[1].tolist().count(0),EpiExp[1].tolist().count(1),EpiExp[1].tolist().count(2)],[EpiExp[2].tolist().count(0),EpiExp[2].tolist().count(1),EpiExp[2].tolist().count(2)],[EpiExp[3].tolist().count(0),EpiExp[3].tolist().count(1),EpiExp[3].tolist().count(2)],[EpiExp[4].tolist().count(0),EpiExp[4].tolist().count(1),EpiExp[4].tolist().count(2)],[EpiExp[5].tolist().count(0),EpiExp[5].tolist().count(1),EpiExp[5].tolist().count(2)]]

sumtable = pandas.DataFrame(sumdata,transpositions,episclasses)
exptable = pandas.DataFrame(expdata,transpositions,episclasses)

print('Reduction by sum: percentages of each class per transposition \n')
print(sumtable/15, '\n')

print('total percentage for each class \n')
print(sumtable.sum(axis=0)/90, '\n')



print('Reduction by exponent: percentages of each class per transposition \n')

print(exptable/15,'\n')
print('total percentage for each class \n')
print(exptable.sum(axis=0)/90, '\n')

