import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as figure
import scipy.spatial as spat

test = np.load('THENEWTHING.npy')
truth = np.load('friday_hack_ABC.npy')


params = np.load('params.npy')

truthmean = truth.mean(axis = 0)
truthstd = truth.std(axis = 0)
eps = 3

count = np.zeros(len(test[:,0,0]))
AllNorms = np.zeros((len(test[:,0,0]),len(test[0,:,0])))

for i in range(len(test[:,0,0])):
    for j  in range(len(test[0,:,0])):
        
        nrm = spat.distance.pdist([test[i,j,:],truthmean])
        if nrm < eps:
            count[i]+=1
        AllNorms[i,j] = nrm  
               


count = count/len(test[0,:,0])


##theres no way this is the best way to do this 
my_xticks = []
for q in range(len(params)):
    my_xticks.append(str(params[q,:]))


plt.xticks(np.arange(0,len(count)),my_xticks, rotation=90)

plt.plot(count)

plt.show()

#[i for i,x in enumerate(count) if x==.2]

