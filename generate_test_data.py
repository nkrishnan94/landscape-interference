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
import itertools

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
		if r <= s:
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


iterations_per_parameter_set = 500
coarseness = 2
count = 0
#this makes the REAL NARA
newthing = np.empty((coarseness**(2**N), iterations_per_parameter_set, 2**N))

a = [range(0,coarseness)]*(2**N)
iterator = tuple(itertools.product(*a))
ref_vec  = np.arange(1,coarseness+1)/coarseness
ref_vec = ref_vec.tolist()
params = []

for thing1 in iterator:
    landscape = [ref_vec[thing1[0]], ref_vec[thing1[1]], ref_vec[thing1[2]], ref_vec[thing1[3]]]
    #save params
    params.append(landscape)


    for j in range(iterations_per_parameter_set):

        population = [0 for i in range(2**N)]
        population[0] = 10
        life_history = [population]

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

            t += 1
            population = copy(new_pop)

            pop = life_history[-1]

        newthing[count,j,:] = pop
    count += 1				
					
np.save('THENEWTHING.npy', newthing)
np.save('params.npy',params)

