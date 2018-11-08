#!/usr/bin/env python

usage='''

chimera --start "Build Structure" --nogui --script "psi_n_phi_angles.py"\n

**NOTE** You must hardcode your mutation info (look inside script)
1) Loads pdb or gro into chimera
2) Changes the psi and phi angles of residues that user specifies
3) Saves new pdb

'''
print usage
#chimera --script "script.py [pdb] [mutation] [residue.chain]"\n

# Import Modules:{{{
import chimera
# }}}

# Methods:{{{

# Rotate Phi and Psi angles():{{{

def rotate_angles(phi,psi,chain,model=0):
    from Midas.midas_text import makeCommand as mc
    model = chimera.openModels.list()[int(model)]
    #IDS = [];
    for res in model.residues:
	ID = res.id
        chimeraID = ':'+str(ID)+'.'+str(chain)
        mc('sel %s'%chimeraID)
        mc('setattr r phi %s sel '%phi)
        mc('setattr r psi %s sel'%psi)


# >>>
rotate_angles(phi='-47.0',psi='-57.8',chain='')
#}}}



