import numpy as np
import matplotlib.pyplot as plt

test = np.load('THENEWTHING.npy')
truth = np.load('friday_hack_ABC.npy')

truthmean = truth.mean(axis = 0)
teststd = test.std(axis = 0)
eps = 50


count = np.zeros(len(test[:,0,0]))
for i in range(len(test[:,0,0])):
    for j  in range(len(test[0,:,0])):
        nrm = np.linalg.norm(test[i,j,:]-truthmean)
        print(nrm)
        if nrm < eps:
            count[i]+=1
                

count = count/len(test[0,:,0])
plt.plot(count)
plt.show()


