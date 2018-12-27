#!/usr/bin bash
gmx editconf -f test.pdb -o newbox.gro -bt cubic -d 1.0

gmx solvate -cp newbox.gro -cs spc216.gro -p new.top -o solv.gro

gmx grompp -f mdp/minimize.mdp -c solv.gro -p new.top -o ions.tpr -maxwarn -1

gmx genion -s ions.tpr -o solv_ions.gro -p new.top -pname NA -nname CL -neutral -conc 0.1

gmx make_ndx -f solv_ions.gro
