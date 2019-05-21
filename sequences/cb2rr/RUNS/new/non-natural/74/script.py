from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm import XmlSerializer

gro =  GromacsGroFile('solvent_ions_equilibrated2.gro')
top = GromacsTopFile('new.top', periodicBoxVectors=gro.getPeriodicBoxVectors())

system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometers, constraints=HBonds)
with open('system.xml','w') as f:
    input=XmlSerializer.serialize(system)
    f.write(input)

integrator = LangevinIntegrator(300*kelvin, 1.0/picoseconds, 2.0*femtoseconds)
with open('integrator.xml','w') as f:
    input=XmlSerializer.serialize(integrator)
    f.write(input)

simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
#simulation.context.setVelocities(velocities)
simulation.context.setVelocitiesToTemperature(300*kelvin)
State = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True, getParameters=True)
with open('state.xml','w') as f:
    input=XmlSerializer.serialize(State)
    f.write(input)
