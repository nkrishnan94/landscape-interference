import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as figure
import scipy.spatial as spat
import statsmodels.tsa.stattools as stattools


test = np.load('THENEWTHING.npy')
truth = np.load('friday_hack_ABC.npy')


params = np.load('params.npy')

truthmean = truth.mean(axis = 0)
truthstd = truth.std(axis = 0)

#this is the tolerance
eps = 32

count = np.zeros(len(test[:,0,0]))
AllNorms = np.zeros((len(test[:,0,0]),len(test[0,:,0])))

for i in range(len(test[:,0,0])):
    for j  in range(len(test[0,:,0])):
        
        #find covariance for each allele
        nrm = abs(np.cov([test[i,j,:,0],truthmean[:,0]]))
        nrm1 = abs(np.cov([test[i,j,:,1],truthmean[:,1]]))
        nrm2 = abs(np.cov([test[i,j,:,2],truthmean[:,2]]))
        nrm3 = abs(np.cov([test[i,j,:,3],truthmean[:,3]]))
        
        #something like explained variation/unexplained variation.
        #each covariance is divided by variance in the allele in the reference set
        AllAlleleCov = [nrm[1,0]/nrm[1,1],nrm1[1,0]/nrm1[1,1],nrm2[1,0]/nrm2[1,1],nrm3[1,0]/nrm3[1,1]]
        
        # rejected if its ABOVE the set tolerance, which is odd, i think.    
        if sum(AllAlleleCov) > eps:
            count[i]+=1
        
        #saving all the summary statistics
        AllNorms[i,j] = sum(AllAlleleCov)



count = count/len(test[0,:,0])


##theres no way this is the best way to do this 
my_xticks = []
for q in range(len(params)):
    my_xticks.append(str(params[q,:]))


plt.xticks(np.arange(0,len(count)),my_xticks, rotation=90)

plt.plot(count)

plt.show()

#[i for i,x in enumerate(count) if x==.2]

