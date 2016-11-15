# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 19:48:35 2016

@author: Nikhil Krishnan
"""

import numpy as np
from copy import copy
from copy import deepcopy
from random import random
from random import choice
import matplotlib.pyplot as plt
import pylab
import os

dir_path = os.path.dirname(os.path.realpath('landscape_inference_parameters.py'))
runfile('simulation.py', dir_path)

#parameter ranges
mut_min=.1
mut_max=.2
mut_step=.1
death_min=.01
death_max=.02
death_step=0.01
drug="AMP"
samples=10

for x in np.arange(mut_min,mut_max,mut_step):
    for y in np.arange(death_min,death_max, death_step):
        simulation(y,x,drug,samples)
        
#function to load and plot a particular result
def viewResult(deathrate,mutrate,drg,samp):
    #load   
    import pandas as pd
    data = pd.read_csv("life_historyavg_"+str(deathrate)+"_"+str(mutrate)+"_"+drg+"_"+str(samp)+".csv")
    data = pd.DataFrame.as_matrix(data)
   
    #plot
    #assume N=4 and T_max =100
    N=4
    T_max=100
    for j in range(2**N):
        plt.plot([data[i][j] for i in range(T_max)],label=1)
        
        
    
    
    
    
        
