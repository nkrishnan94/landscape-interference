import numpy as np
import matplotlib.pyplot as plt
import matplotlib.figure as figure
import scipy.spatial as spat
from mpl_toolkits.axes_grid1 import AxesGrid





test = np.load('THENEWTHING.npy')
truth = np.load('friday_hack_ABC.npy')


params = np.load('params.npy')

truthmean = truth.mean(axis = 0)
truthstd = truth.std(axis = 0)

#this is the tolerance
eps = 75000
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
        if nrm[1,0]> eps:
            count[i]+=1
            
        #saving all the summary statistics
        AllNorms[i,j] = sum(AllAlleleCov)



count = count/len(test[0,:,0])


#theres no way this is the best way to do this 

my_xticks = []
for q in range(len(params)):
    my_xticks.append(str(params[q,:]))


plt.xticks(np.arange(0,len(count)),my_xticks, rotation=90)

plt.plot(count)

plt.show()

[i for i,x in enumerate(count) if x==.2]


trace = count.reshape((len(count)**.25),(len(count)**.25),(len(count)**.25),(len(count)**.25))
#f, axarr = plt.subplots(int(len(count)**.25), int(len(count)**.25))
#for q in range(len(params)):
#for k in range(int(len(count)**.25)):
#    for l in range(int(len(count)**.25)):
#        plt.imshow(trace[k,l], cmap='hot', interpolation='nearest')
        
#fig = plt.figure()
#vals = count.reshape((len(count)**.5),(len(count)**.25),(len(count)**.25))
#grid = AxesGrid(fig,111,nrows_ncols=(2, 2),axes_pad=0.5,label_mode = 'all')


ti = 0;
it= len(count)**.25;
#for val, ax in zip(vals,grid):
#
#    im = ax.imshow(val)
#    ax.set_xticklabels(['',.5,'',1,''])
#    ax.set_yticklabels(['',.5,'',1,''])
#    if ti < int(len(count)**.25):
#        ax.set_title('allele 2 = %s'%(params[ti,3]))
#    if it % len(count)**.25 == 0:
#         ax.set_ylabel('allele 1 = %s'%(params[ti,3]))
#    ti =ti +1
#    it= it +1

    


 