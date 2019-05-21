import sys, os
import numpy as np
seq=['A8','A8_M','A13','A13_M','A18','A18_M']
nRUN = [3,4,3,3,3,3]
proj=np.arange(14140,14146,1)
for i in range(len(proj)):
	for j in range(nRUN[i]):
		os.chdir('for_vav4/%d/RUN%d'%(proj[i],j))
		os.system('python script.py')
		os.chdir('../../../')





