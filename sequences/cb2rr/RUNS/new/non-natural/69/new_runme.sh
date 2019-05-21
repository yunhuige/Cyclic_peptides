#!/usr/bin bash
gmx editconf -f ligand_GMX_invert.pdb -o newbox.gro -bt cubic -c -box 4.06871 4.06871 4.06871

gmx solvate -cp newbox.gro -cs spc216.gro -p new.top -o solv.gro

gmx grompp -f mdp/minimize.mdp -c solv.gro -p new.top -o ions.tpr -maxwarn -1

gmx genion -s ions.tpr -o solv_ions.gro -p new.top -pname NA -nname CL -neutral -conc 0.1

gmx make_ndx -f solv_ions.gro
