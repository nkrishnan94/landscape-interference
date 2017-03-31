# nothing here yet# Daniel Nichol and Jacob G Scott
# CA to simulate bacteria experiments.
# 15/4/2015
import numpy as np
#import scipy as sp
from copy import copy
from copy import deepcopy
from random import randint
from random import random
from random import choice


#Simulation constants.

p_divide = 1.0 #Probability to divide per time-step (time-step = 30mins, typical doubling time.)
N = 2 #Genotype length
cc = 10**3 #carrying capacity
death_prob = 0.15 #probability of death per timestep.
mutation_prob = 10**-2
T_max = 50
drug_conc = 0

# reproduction probability
def get_p_divide(MIC, drug_conc):
	prob = (1. - (drug_conc/MIC)) * p_divide
	return prob

def choose_random(pop):
	r = randint(1,sum(pop))
	s = 0
	gt = 0
	for i in range(len(pop)):
		s += pop[i]
		if r<=s:
			gt = i
			break
	return gt

	#============================== Helper Functions ==============================#

# Takes a genotype and converts it to an integer use for indexing lists
def convertGenotypeToInt(genotype):
		out = 0
		for bit in genotype:
			out = (out << 1) | bit
		return out

# Converts an integer to a genotype by taking the binary value and padding to the left by 0s		
def unindex(anInt, pad):
	offset = 2**pad
	return [int(x) for x in bin(offset+anInt)[3:]]	

# Function which returns all genotypes at Hamming distance 1 from a specified genotype
def oneStepNeighbours(genotype):
	neighbours = []
	for x in range(0, len(genotype)):
		temp = deepcopy(genotype)
		temp[x] = (genotype[x]+1) %2 #There is some inefficiency here.
		neighbours.append(temp)
	return neighbours	

	#============================== Landscapes ==============================#

NARA = [0.5, 0.5, 0.5, 1.0]

#AMP = [1.851, 2.082, 1.948, 2.434, 2.024, 2.198, 2.033, 0.034, 1.57, 2.165, 0.051, 0.083, 2.186, 2.322, 0.088, 2.821]
# AM  = [1.778, 1.782, 2.042, 1.752, 1.448, 1.544, 1.184, 0.063, 1.72, 2.008, 1.799, 2.005, 1.557, 2.247, 1.768, 2.047]
# CEC	= [2.258, 1.996, 2.151, 2.648, 2.396, 1.846, 2.23, 0.214, 0.234, 0.172, 2.242, 0.093, 2.15, 0.095, 2.64, 0.516]
# CTX = [0.16, 0.085, 1.936, 2.348, 1.653, 0.138, 2.295, 2.269, 0.185, 0.14, 1.969, 0.203, 0.225, 0.092, 0.119, 2.412]
# ZOX = [0.993, 0.805, 2.069, 2.683, 1.698, 2.01, 2.138, 2.688, 1.106, 1.171, 1.894, 0.681, 1.116, 1.105, 1.103, 2.591]
# CXM = [1.748, 1.7, 2.07, 1.938, 2.94, 2.173, 2.918, 3.272, 0.423, 1.578, 1.911, 2.754, 2.024, 1.678, 1.591, 2.923]
# CRO = [1.092, 0.287, 2.554, 3.042, 2.88, 0.656, 2.732, 0.436, 0.83, 0.54, 3.173, 1.153, 1.407, 0.751, 2.74, 3.227]
# AMC = [1.435, 1.573, 1.061, 1.457, 1.672, 1.625, 0.073, 0.068, 1.417, 1.351, 1.538, 1.59, 1.377, 1.914, 1.307, 1.728]
# CAZ = [2.134, 2.656, 2.618, 2.688, 2.042, 2.756, 2.924, 0.251, 0.288, 0.576, 1.604, 1.378, 2.63, 2.677, 2.893, 2.563]
# CTT = [2.125, 1.922, 2.804, 0.588, 3.291, 2.888, 3.082, 3.508, 3.238, 2.966, 2.883, 0.89, 0.546, 3.181, 3.193, 2.543]
# SAM = [1.879, 2.533, 0.133, 0.094, 2.456, 2.437, 0.083, 0.094, 2.198, 2.57, 2.308, 2.886, 2.504, 3.002, 2.528, 3.453]
# CPR = [1.743, 1.662, 1.763, 1.785, 2.018, 2.05, 2.042, 0.218, 1.553, 0.256, 0.165, 0.221, 0.223, 0.239, 1.811, 0.288]
# CPD = [0.595, 0.245, 2.604, 3.043, 1.761, 1.471, 2.91, 3.096, 0.432, 0.388, 2.651, 1.103, 0.638, 0.986, 0.963, 3.268]
# TZP = [2.679, 2.906, 2.427, 0.141, 3.038, 3.309, 2.528, 0.143, 2.709, 2.5, 0.172, 0.093, 2.453, 2.739, 0.609, 0.171]
# FEP = [2.59, 2.572, 2.393, 2.832, 2.44, 2.808, 2.652, 0.611, 2.067, 2.446, 2.957, 2.633, 2.735, 2.863, 2.796, 3.203]

landscape = NARA

# landscape = map(np.exp, landscape)

def update(pop, g_index, drug_conc):
	r1,r2,r3 = random(), random(), random()
	if pop[g_index] == 0:
		print("WTF IS GOING ON!?")
	delta_pop = [0 for i in range(16)]
	new_gi = g_index
	if r1 < death_prob:
		delta_pop[g_index] = -1
	elif sum(pop) < cc and r2<get_p_divide(landscape[g_index], drug_conc):
		if r3 < mutation_prob:
			osn = oneStepNeighbours(unindex(g_index,N))
			new_g = choice(osn)
			new_gi = convertGenotypeToInt(new_g)
		delta_pop[new_gi]+=1

	for k in range(len(pop)):
		pop[k] = pop[k] + delta_pop[k]
	return pop


# def simpson_index(pop):
# 	tot = float(sum(pop))
# 	sq_props = map(lambda x : (x/tot)**2, pop)
# 	si = sum(sq_props)
# 	return 1-si

# def shannon_entropy(pop):
# 	tot = float(sum(pop))
# 	shi = 0
# 	for p in pop:
# 		if p!=0:
# 			shi+=(p/tot)*np.log(p/tot)
# 	return -shi

# def get_allele_frequencies(population):
# 	afs = [0.,0.,0.,0.]
# 	for i in range(2**N):
# 		b = unindex(i,N)
# 		for j in range(len(b)):
# 			if b[j]==1:
# 				afs[j]+=population[i]
# 	return np.array(afs)/sum(population)


#this makes the REAL NARA
thing = []

for j in range(1000):
	#Recording information
	population = [0 for i in range(2**N)]
	population[0]=10
	life_history = [population]
	# afs_history = [get_allele_frequencies(population)]
	# div_history = [simpson_index(population)]
	# ent_history = [shannon_entropy(population)]

	t = 0
	while t < T_max:
		# print t
		#plt.clf()
		new_pop = copy(population)
		count_list = copy(population)
		while sum(count_list)!=0:
			g_index = choose_random(count_list)
			count_list[g_index]-=1
			new_pop = update(new_pop, g_index, drug_conc)
		life_history.append(population)
		# div_history.append(simpson_index(population))
		# ent_history.append(shannon_entropy(population))
		# afs_history.append(get_allele_frequencies(population))
		t+=1
		population = copy(new_pop)

		pop = life_history[-1]
	thing.append(pop)
	
np.save('friday_hack_ABC.npy', thing)

