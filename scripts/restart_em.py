from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import pandas as pd
import openmmtools

# Set parameters
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

# Load checkpoint
simulation.loadCheckpoint(output_prefix + '.chk')

# Set up reporters for state data, checkpoint file, and trajectory 
simulation.reporters.append(StateDataReporter(output_prefix + '.restarted.csv', 12500, time=True, step=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True, volume=True, temperature=True))
simulation.reporters.append(CheckpointReporter(output_prefix + '.restarted.chk', 12500))
simulation.reporters.append(DCDReporter(output_prefix +  '.restarted.dcd', 12500, enforcePeriodicBox=False))

# Check last time step
df_data = pd.read_csv(output_prefix + '.csv')
last_step = df_data.iloc[-1]['#"Step"']

# Start simulation
simulation.step(12500000 - last_step)

