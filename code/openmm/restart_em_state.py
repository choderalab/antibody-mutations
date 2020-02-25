from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import pandas as pd
import openmmtools
import time
import progressbar
from simtk.openmm import XmlSerializer
import math

# Set parameters
timestep = 4*femtosecond
input_file = sys.argv[1]
output_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_output/" + os.path.basename(input_file)[:-18] + '50ns'

# Set up forcefield
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Read in energy minimized PDB
pdb = PDBFile(input_file)

# Deserialize system file and load system
with open(output_prefix + '.xml', 'r') as f:
    system = XmlSerializer.deserialize(f.read())

# Set up integrator
integrator = openmmtools.integrators.LangevinIntegrator(300*kelvin, 1/picosecond, 4*femtosecond)

# Set up platform
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Set up simulation 
simulation = Simulation(pdb.topology, system, integrator, platform)

# Load state and set box vectors, positions, and velocities
with open(output_prefix + '.state.xml', 'r') as infile:
    state_xml = infile.read()
state = XmlSerializer.deserialize(state_xml)
simulation.context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
simulation.context.setPositions(state.getPositions())
simulation.context.setVelocities(state.getVelocities())
simulation.context.setTime(state.getTime())

# Set up reporters for state data, checkpoint file, and trajectory 
simulation.reporters.append(StateDataReporter(output_prefix + '.restarted.csv', 12500, time=True, step=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True, volume=True, temperature=True))
simulation.reporters.append(DCDReporter(output_prefix +  '.restarted.dcd', 12500, enforcePeriodicBox=False))

# Check last time step
df_data = pd.read_csv(output_prefix + '.csv')
last_step = df_data.iloc[-1]['#"Step"']
nsteps = 12500 # 50 ps
niterations = math.ceil((12500000 - last_step)/12500)

# Equilibrate
print('Equilibrating...')
initial_time = time.time()
for iteration in progressbar.progressbar(range(niterations)):
    simulation.step(nsteps)
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(output_prefix + '.state.restarted.xml', 'w') as outfile:
        state_xml = XmlSerializer.serialize(state)
        outfile.write(state_xml)
elapsed_time = (time.time() - initial_time) * seconds
simulation_time = niterations * nsteps * timestep
print('    Equilibration took %.3f s for %.3f ns (%8.3f ns/day)' % (elapsed_time / seconds, simulation_time / nanoseconds, simulation_time / elapsed_time * day / nanoseconds))

# Save equilibrated pdb
print("Writing equilbrated model to ", output_prefix + ".equilibrated.pdb")
positions_equilibrated = simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions()
PDBFile.writeFile(modeller.topology, positions_equilibrated, open(output_prefix + ".equilibrated.pdb", 'w'), keepIds=True)
print('  final   : %8.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy()/kilocalories_per_mole))

