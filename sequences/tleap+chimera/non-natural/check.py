import sys, os

for i in range(63):
    print os.path.isfile('%d/solv_ions.gro'%i)
