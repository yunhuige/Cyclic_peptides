#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -N ligand_test
#PBS -q gpu
#PBS -l nodes=1:ppn=1
#PBS -o ligand_test 
#PBS 

cd $PBS_O_WORKDIR

module load cuda

python test_openmm.py
