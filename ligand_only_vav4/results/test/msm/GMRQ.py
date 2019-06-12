import os,sys
import pandas as pd
import numpy as np
import pickle
from msmbuilder.cluster import KCenters, KMeans
from sklearn.cross_validation import KFold
from msmbuilder.msm import MarkovStateModel
import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt


tica_data = np.load('tica.npy')

cv = KFold(len(tica_data),n_folds=5)
results = []


#n_clusters = [100,200,400,600,800,1000]
n_clusters = [5,10,25,50,75,100,150, 200]
#n_clusters = [10,100,200,400,600]
#n_clusters = [2500,3000,3500]
#n_clusters=[5,10]
lagtime = 400

for n in n_clusters:
    kcenters = KCenters(n_clusters=n)
    print "Clustering data to %d clusters..."%n
    for fold,(train_index,test_index) in enumerate(cv):
#        reduced_train_index = train_index[::100]
#	print 'train_index:', train_index,'test_index:', test_index
        kcenters.fit(tica_data[train_index])
        assignments_train = kcenters.transform(tica_data[train_index])
        assignments_test = kcenters.transform(tica_data[test_index])
        msm = MarkovStateModel(n_timescales=10,lag_time=lagtime)
        msm.fit(assignments_train)
        train_score = msm.score(assignments_train)
        test_score = msm.score(assignments_test)
        if test_score > 0.0:
            results.append({'train_score':train_score,
                            'test_score':test_score,
                            'n_states':n,
                            'fold':fold})

#        results.append({'train_score':train_score,
#                        'test_score':test_score,
#                        'n_states':n,
#                        'fold':fold})
#sys.exit()
results = pd.DataFrame(results)
print results
np.save('GMRQ_test.npy',results)
#print results.head()

plt.figure(figsize=(8,6))
plt.plot(results['n_states'], results['train_score'], c='b', marker='.', ls='')
plt.plot(results['n_states'], results['test_score'], c='r', marker='.', ls='')

mean_over_folds = results.groupby('n_states').aggregate(np.mean)
plt.plot(mean_over_folds.index, mean_over_folds['test_score'], c='r', marker='.', ls='-', label='Mean test')
plt.plot(mean_over_folds.index, mean_over_folds['train_score'], c='b', marker='.', ls='-', label='Mean train')
plt.semilogx()
plt.ylabel('Generalized Matrix Rayleigh Quotient (Score)')
plt.xlabel('Number of states')

best_n_states = np.argmax(mean_over_folds['test_score'])
best_test_score = mean_over_folds.ix[best_n_states]['test_score']
plt.plot(best_n_states, best_test_score, marker='*', ms=20, c='w', label='n_states=%d' % best_n_states)

plt.legend(loc='best', numpoints=1)
plt.savefig('test_GMRQ.png')
plt.show()

