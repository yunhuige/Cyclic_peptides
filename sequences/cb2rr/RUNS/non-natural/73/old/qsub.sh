#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -N ligand
#PBS -q normal
#PBS -l nodes=1:ppn=20
#PBS -o ligand 
#PBS 

cd $PBS_O_WORKDIR

module load gromacs/5.1.0
module load openmpi

gmx grompp -f mdp/minimize.mdp -c solv_ions.gro -p new.top -n index.ndx -o minimize2.tpr -maxwarn -1
gmx mdrun -v -nt 1 -s minimize2.tpr -c solvent_ions_minimized.gro
gmx grompp -f mdp/equil.mdp -c solvent_ions_minimized.gro -p new.top -n index.ndx -o equil.tpr -maxwarn -1
gmx mdrun -v -nt 2 -s equil.tpr -c solvent_ions_equilibrated.gro
gmx grompp -f mdp/equil_NVT.mdp -c solvent_ions_equilibrated.gro -p new.top -n index.ndx -o equil2.tpr -maxwarn -1
gmx mdrun -v -nt 2 -s equil2.tpr -c solvent_ions_equilibrated2.gro
