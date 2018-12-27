import sys, os
for i in range(63,74):
	os.system('mkdir -p for_vav4/%d'%i)
	os.chdir('%d/'%i)
	os.system('cp -r amber99sbnmr1-ildn.ff solvent_ions_equilibrated2.gro new.top ../for_vav4/%d'%i)
	os.chdir('../')
	os.system('cp script.py for_vav4/%d'%i)
