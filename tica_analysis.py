import sys, os
import numpy as np
from itertools import combinations
import matplotlib
matplotlib.use('Agg') # must be before pyplot as plt
from matplotlib import pyplot as plt
fontfamily={'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


#NOTE: Does n=5 ???
def triangle_subplots(data,n=5):
    print("Plotting triangle subplots for tICA...")
    loc = [1,5,9,13,6,10,14,11,15,16]
    combos = list(combinations(range(1,int(n)+1),2))
    fig = plt.figure(figsize=(12,10))
    for i in range(len(loc)):
        ax = fig.add_subplot(4,4,int('%i'%(loc[i])))

        ax.hist2d(data[:,int(combos[i][0]-1)], data[:,int(combos[i][1]-1)],
                bins=200, norm=LogNorm(), cmap="jet")

        if loc[i] in [13,14,15,16]:
            ax.set_xlabel(r'$tIC_{%s}$'%(str(combos[i][0])), fontsize=20)
        if loc[i] in [1,5,9,13]:
            ax.set_ylabel(r'$tIC_{%s}$'%(str(combos[i][1])), fontsize=20)

    ticks = [ax.xaxis.get_minor_ticks(),ax.yaxis.get_minor_ticks(),
             ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()]
    for k in range(0,len(ticks)):
        if k == 0 or 1:
            for tick in ticks[k]:
                tick.label.set_fontsize(16)
        else:
            for tick in ticks[k]:
                tick.label.set_fontsize(18)
    plt.tight_layout()
    plt.savefig('tICA_triangle.pdf')


def plot_tICA(data, tICs=[0,1], output='tIC_1_and_2'):
    print('Plotting tICA...')
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111)
    # Sometimes the sign of the tICA components need to be flipped
    ax.hist2d(data[:,tICs[0]], data[:,tICs[1]], bins=200, norm=LogNorm(), cmap="jet")
    ax.set_xlabel(r'$tIC_{1}$', fontsize=20)
    ax.set_ylabel(r'$tIC_{2}$', fontsize=20)
    ticks = [ax.xaxis.get_minor_ticks(),ax.yaxis.get_minor_ticks(),
             ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()]
    for k in range(0,len(ticks)):
        if k == 0 or 1:
            for tick in ticks[k]:
                tick.label.set_fontsize(16)
        else:
            for tick in ticks[k]:
                tick.label.set_fontsize(18)
    plt.tight_layout()
    plt.savefig('%s.pdf'%output)



def plot_traces(data, traj_path, tICs=[0,1], output='tIC_1_and_2_traces'):
    print("Plotting traces on tICA...")
    tica=np.load('%s'%(traj_path))
    for k in range(len(tica)):
        print(k)
        t1, t2 = [], []
        for l in range(len(tica[k])):
            t1.append(tica[k][l][0])
            t2.append(tica[k][l][1])
        x=np.arange(0,len(t1)/10.0,0.1)
        plt.subplot(2,2,1)
        plt.hist2d(data[:,tICs[0]], -data[:,tICs[1]],
                bins=300, norm=LogNorm(), cmap="jet")
        plt.xlabel('tIC1',fontsize=8)
        plt.ylabel('tIC2',fontsize=8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.plot(t1,t2,color='black',ms=1,linewidth=0.5)
        plt.plot(t1[0],t2[0],'*',color='magenta',ms=10.0)
        plt.plot(t1[-1],t2[-1],'*',color='green',ms=10.0)
        plt.xlim(-1,8)
        plt.ylim(-15,5)
        plt.title('Pathway traj%d'%k,fontsize=10)
        plt.subplot(2,2,3)
        plt.plot(t1,x)
        plt.xlabel('tIC1',fontsize=8)
        plt.ylabel('time (ns)',fontsize=8)
        plt.xlim(-1,8)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.subplot(2,2,2)
        plt.plot(x,t2)
        plt.xlabel('time (ns)',fontsize=8)
        plt.ylabel('tIC2',fontsize=8)
        plt.ylim(-15,5)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
        plt.tight_layout()
        plt.savefig('%s.pdf'%output)



for s in seq:
    os.chdir('%s'%s)
    a=np.load('tica.npy')
    b=[]
    c=[]
    d=[]
    e=[]
    for i in range(len(a)):
        for j in range(len(a[i])):
            b.append(a[i][j][0])
            c.append(a[i][j][1])
            d.append(a[i][j][2])
            e.append(a[i][j][3])
    np.save('tica1.npy',b)
    np.save('tica2.npy',c)
    np.save('tica3.npy',d)
    np.save('tica4.npy',e)

    tica_data = np.array([b,c,d,e])
    triangle_subplots(tica_data)
    plot_tICA(tica_data)
    traj_path = "../numpyfile for trajectory to be overlaied"
    plot_traces(tica_data, traj_path)


