from __future__ import print_function
import os,sys
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

class Tee(object):
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)


gro = GromacsGroFile('solvent_ions_equilibrated2.gro')
#gro = GromacsGroFile('solv_ions.gro')
top = GromacsTopFile('new.top', periodicBoxVectors=gro.getPeriodicBoxVectors())

system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

simulation = Simulation(top.topology, system, integrator)
#simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)

platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaPrecision': 'mixed','CudaDeviceIndex': gpu_id}
properties = {'CudaPrecision': 'mixed'}

simulation.reporters.append(DCDReporter('traj.dcd', 50000))
simulation.reporters.append(StateDataReporter(Tee('tarj.log', 'w'), 50000, step=True, potentialEnergy=True, temperature=True, speed=True, separator='\t'))
simulation.step(2500000)
#simulation.runForClockTime(22,checkpointFile="checkpoint.chk")
simulation.saveCheckpoint("checkpoint.chk")
sys.exit()

print("\nInitial system energy")
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy()
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

a=simulation.context.getState(getPositions=True)

#starting_pdb=sys.argv[1]
#sequence_number=int(sys.argv[2])

#pdb = app.PDBFile(starting_pdb)
#forcefield = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')

#system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
#    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True)
#integrator = mm.LangevinIntegrator(300*unit.kelvin, 91.0/unit.picoseconds, 
#    2.0*unit.femtoseconds)
#integrator.setConstraintTolerance(0.00001)

platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaPrecision': 'mixed','CudaDeviceIndex': gpu_id}
properties = {'CudaPrecision': 'mixed'}
integrator2 = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

simulation2 = Simulation(top.topology, system, integrator2)
#simulation = Simulation(top.topology, system, integrator)
simulation2.context.setState(a)

#print('Equilibrating...')
#simulation.step(1000)

#print("\nInitial system energy")
#print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
#simulation.minimizeEnergy(maxIterations=200000)
#print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
#initial_positions = simulation.context.getState(getPositions=True).getPositions()
#interval = int(1*picosecond / (2*femtosecond))
#n_steps = int(100*nanosecond / (2*femtosecond))
#interval=50000
#n_steps=1000000

#print('Equilibrating...')
#simulation2.step(1000)

#initial_positions = simulation2.context.getState(getPositions=True).getPositions()

simulation2.reporters.append(DCDReporter('traj.dcd', 1000))
simulation2.reporters.append(StateDataReporter(Tee('tarj.log', 'w'), 1000, step=True, potentialEnergy=True, temperature=True, speed=True, separator='\t'))
#simulation.step(n_steps)
simulation2.runForClockTime(22,checkpointFile="checkpoint.chk")
simulation2.saveCheckpoint("checkpoint.chk")
