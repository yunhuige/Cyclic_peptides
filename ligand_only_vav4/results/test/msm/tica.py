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

distances=[]
for k in range(100):
	print k
        b=np.load('distance%d.npy'%k)
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


