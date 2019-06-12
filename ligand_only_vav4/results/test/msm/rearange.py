import numpy as np
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

