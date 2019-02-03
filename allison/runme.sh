#!/usr/bin bash
gmx editconf -f ligand_GMX_invert.pdb -o newbox.gro -bt cubic -d 1.0

gmx solvate -cp newbox.gro -p ligand_GMX.top -o solv.gro

gmx grompp -f mdp/shortmin.mdp -c solv.gro -p ligand_GMX.top -o ions.tpr

gmx genion -s ions.tpr -o solv_ions.gro -p ligand_GMX.top -pname NA -nname CL -neutral -conc 0.1

gmx make_ndx -f solv_ions.gro


