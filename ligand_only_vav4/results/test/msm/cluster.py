from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.dataset import dataset
import numpy as np
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
import mdtraj as md
import os,sys
import msmbuilder.utils
from msmbuilder.utils import load
from msmbuilder.cluster import KMeans
from msmbuilder.cluster import KCenters
from msmbuilder.cluster import KMedoids
from msmbuilder.cluster import MiniBatchKMeans
from msmbuilder.msm import implied_timescales
from msmbuilder.msm import ContinuousTimeMSM, MarkovStateModel
from itertools import combinations
from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.decomposition import tICA
from sklearn.pipeline import Pipeline
from msmbuilder.example_datasets import fetch_met_enkephalin
from matplotlib import pyplot as plt
from sklearn.externals import joblib
#Featurization
#trajs = md.load('traj0.xtc',top='conf.gro')
#Ind = trajs.topology.select("all")
#trajs1=trajs.atom_slice(Ind)

#distances=np.load('proj8609.npy')
#tica_transformed = tICA(lag_time=10,n_components=4).fit_transform(distances)
#np.save('tica.npy',tica_transformed)
#sys.exit()
tica_transformed=np.load('tica.npy')
#assignments=KCenters(n_clusters=25).fit_predict(tica_transformed)
assignments=KCenters(n_clusters=25).fit(tica_transformed)
centers=assignments.cluster_centers_
np.save('cluster_centers.npy',centers)
ass=assignments.transform(tica_transformed)
#sys.exit()
#for i in range(100):
#        traj2[i].save_pdb('Gens%d.pdb'%i)
np.save('Assignments.npy',ass)
#sys.exit()
lagtimes = np.array([1,5,10,15,30,100,150,200,300,600,1000,2000])
msmts = []
for lagtime in lagtimes:
	print "\tLagtime: %d"%lagtime
    	msm = MarkovStateModel(lag_time=lagtime).fit(ass)
    	msmts.append(msm.timescales_)
np.save('test_msmts.npy', msmts)


#sys.exit()


ns_per_frame = 0.5
import matplotlib.pyplot as plt
plt.figure()
for j in range(9):
	print j
	plt.plot(ns_per_frame*lagtimes,[ns_per_frame*msmts[i][j] for i in range(len(lagtimes))])
#plt.plot(range(len(lagtimes)),[msmts[i][0] for i in range(len(lagtimes))])
#plt.plot(range(len(lagtimes)),[msmts[i][1] for i in range(len(lagtimes))])
#plt.plot(range(len(lagtimes)),[msmts[i][2] for i in range(len(lagtimes))])
#plt.plot(range(len(lagtimes)),[msmts[i][3] for i in range(len(lagtimes))])
#plt.plot(range(len(lagtimes)),[msmts[i][4] for i in range(len(lagtimes))])
plt.xlabel('lagtime (ns)')
plt.ylabel('implied timescale (ns)')
plt.yscale('log')
plt.savefig('lagtime_50_transform.pdf')
plt.show()


#print "Preparation done, now begin clustering..."
#Cluster
#kcenters=KCenters(n_clusters=250,metric='rmsd').fit(trajs1)
#traj2=kcenters.cluster_centers_
#sequences=kcenters.fit_transform(trajs)
#print(kcenters.summarize())
#assignments=kcenters.predict(trajs1)
#np.save('Assignments.npy',assignments)
#for i in range(250):
#        traj2[i].save_pdb('Gens%d.pdb'%i)

