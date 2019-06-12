import sys, os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

t1=np.load('tica1.npy')
t2=np.load('tica2.npy')
t3=np.load('tica3.npy')
t4=np.load('tica4.npy')

plt.figure(figsize=(14,14))
plt.subplot(4,4,1)
plt.hist2d(t1,t1,bins=300,norm=LogNorm())
#plt.xlabel('tIC1')
plt.ylabel('tIC1')

plt.subplot(4,4,5)
plt.hist2d(t1,t2,bins=300,norm=LogNorm())
#plt.xlabel('tIC1')
plt.ylabel('tIC2')

plt.subplot(4,4,9)
plt.hist2d(t1,t3,bins=300,norm=LogNorm())
#plt.xlabel('tIC1')
plt.ylabel('tIC3')

plt.subplot(4,4,13)
plt.hist2d(t1,t4,bins=300,norm=LogNorm())
plt.xlabel('tIC1')
plt.ylabel('tIC4')

plt.subplot(4,4,6)
plt.hist2d(t2,t2,bins=300,norm=LogNorm())
#plt.xlabel('tIC2')
#plt.ylabel('tIC2')

plt.subplot(4,4,10)
plt.hist2d(t2,t3,bins=300,norm=LogNorm())
#plt.xlabel('tIC2')
#plt.ylabel('tIC3')

plt.subplot(4,4,14)
plt.hist2d(t2,t4,bins=300,norm=LogNorm())
plt.xlabel('tIC2')
#plt.ylabel('tIC4')

plt.subplot(4,4,11)
plt.hist2d(t3,t3,bins=300,norm=LogNorm())
#plt.xlabel('tIC3')
#plt.ylabel('tIC3')

plt.subplot(4,4,15)
plt.hist2d(t3,t4,bins=300,norm=LogNorm())
plt.xlabel('tIC3')
#plt.ylabel('tIC4')

plt.subplot(4,4,16)
plt.hist2d(t4,t4,bins=300,norm=LogNorm())
plt.xlabel('tIC4')
#plt.ylabel('tIC4')
plt.tight_layout()
plt.savefig('tica_total.pdf')

plt.show()
