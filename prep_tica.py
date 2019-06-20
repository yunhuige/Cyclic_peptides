#!/usr/bin/env python
import sys, os
import numpy as np
import mdtraj as md
from msmbuilder.featurizer import DihedralFeaturizer
from msmbuilder.dataset import dataset
import numpy as np
from matplotlib import pyplot as plt
import mdtraj as md
import os,sys, glob
import msmbuilder.utils
import pickle
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


for s in seq:
    os.chdir('%s'%s)
    for i in range(100):
        print(i)
        os.chdir('CLONE%d'%i)
        os.system('gmx trjcat -f traj* -o total.xtc')
        os.chdir('../')


    a=md.load('test.gro')
    ind1 = a.topology.select('name == "CA"')
    ind2 = a.topology.select('name  == "CB"')

    ind = list(ind1)+list(ind2)

    pairs = []
    for i in range(0,len(ind)-1):
        for j in range(i+1,len(ind)):
            pairs.append((ind[i],ind[j]))
    np.save('indices_CA_or_CB.npy',pairs)
    os.system('mkdir distance_CA_or_CB')

    ind=np.load('indices_CA_or_CB.npy')

    for i in range(100):
        print(i)
        t=md.load('CLONE%d/total.xtc'%i, top='test.gro')
        distance=md.compute_distances(t,ind)
        np.save('distance_CA_or_CB/distance%d.npy'%i,distance)

    distances=[]
    for k in range(100):
        print(k)
        b=np.load('distance_CA_or_CB/distance%d.npy'%k)
        distances.append(b)
    tica_model = tICA(lag_time=50,n_components=4)
    tica_fit = tica_model.fit(distances)

    output =open('tica_model.npy','wb')
    pickle.dump(tica_model,output)
    output.close()

    tica_transformed = tica_model.transform(distances)

    output = open('tica_model.eigenvectors.npy', 'wb')
    pickle.dump(tica_model.eigenvectors_, output)
    output.close()

    np.save('tica.npy',tica_transformed)
    os.chdir('../')


