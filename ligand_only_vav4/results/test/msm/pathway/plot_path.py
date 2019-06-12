import sys, os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

a=np.load('../tica1.npy')
b=np.load('../tica2.npy')
tica=np.load('../tica.npy')
for k in range(len(tica)):
#for k in range(1):
        print k
        t1=[]
        t2=[]
        for l in range(len(tica[k])):
                t1.append(tica[k][l][0])
                t2.append(tica[k][l][1])
	x=np.arange(0,len(t1)/10.0,0.1)
#        plt.figure()
	plt.subplot(2,2,1)
        plt.hist2d(a,b,bins=300,norm=LogNorm())
#        plt.colorbar()
        plt.xlabel('tIC1',fontsize=8)
        plt.ylabel('tIC2',fontsize=8)
	plt.xticks(fontsize=8)
	plt.yticks(fontsize=8)
        plt.plot(t1,t2,color='black',ms=1,linewidth=0.5)
#       for i in range(100):
#               plt.plot(d[i],e[i],"o",color='red')
#               plt.annotate('%s'%i,xy=(d[i],e[i]))

        plt.plot(t1[0],t2[0],'*',color='magenta',ms=10.0)
        plt.plot(t1[-1],t2[-1],'*',color='green',ms=10.0)
        plt.xlim(-2,2)
        plt.ylim(-1.5,3)
        plt.title('Pathway traj%d'%k,fontsize=10)
	plt.subplot(2,2,3)
	plt.plot(t1,x)
	plt.xlabel('tIC1',fontsize=8)
	plt.ylabel('time (ns)',fontsize=8)
	plt.xlim(-2,2)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
	plt.subplot(2,2,2)
	plt.plot(x,t2)
	plt.xlabel('time (ns)',fontsize=8)
	plt.ylabel('tIC2',fontsize=8)
	plt.ylim(-1.5,3)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
	plt.tight_layout()
#	plt.savefig('junk.png')
        plt.savefig('path_%d.png'%k)
#        plt.show()
        plt.close()
#       plt.show()

