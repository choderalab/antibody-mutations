from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import openmmtools
import time
import progressbar
from simtk.openmm import XmlSerializer

# Set global parameters
timestep = 4*femtosecond
nsteps = 12500 # 50 ps
niterations = 1000 # 50 ns

# Set file names
input_file = sys.argv[1]
output_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_output/" + os.path.basename(input_file)[:-2] + '50ns'

# Set up forcefield
print("Loading forcefield...")
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Set up topology
print("Loading ", input_file)
pdb = PDBFile(input_file)
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.4)

# Add solvent
print("Adding solvent...")
modeller.addSolvent(forcefield, ionicStrength=150*millimolar, padding=1.0*nanometers) # default: tip3p

# Save solvated PDB
print("Writing the solvated model to ", output_prefix + ".solvated.pdb")
PDBFile.writeFile(modeller.topology, modeller.positions, open(output_prefix + ".solvated.pdb", 'w'), keepIds=True)

# Set up system
print('Creating OpenMM System...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds, hydrogenMass=4*amu)

# Add barostat
print('Adding barostat...')
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

# Set up integrator
integrator = openmmtools.integrators.LangevinIntegrator(300*kelvin, 1/picosecond, timestep)

# Serialize and save the system to an xml file
print("Seralizing the system to ", output_prefix + ".xml")
with open(output_prefix + '.xml', 'w') as f:
	f.write(XmlSerializer.serialize(system))

# Set up platform
print("Setting up the platform...")
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Set up simulation 
print("Setting up the simulation...")
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# Minimize the energy
print("Minimizing the energy...")
print('  initial : %8.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy()/kilocalories_per_mole))
simulation.minimizeEnergy()
print('  final   : %8.3f kcal/mol' % (simulation.context.getState(getEnergy=True).getPotentialEnergy()/kilocalories_per_mole))

# Save energy minimized PDB
print("Writing the minimized model to ", output_prefix + ".minimized.pdb")
positions_minimized = simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions()
PDBFile.writeFile(modeller.topology, positions_minimized, open(output_prefix + ".minimized.pdb", 'w'),  keepIds=True)

# Set up reporters for state data, checkpoint file, and trajectory 
simulation.reporters.append(StateDataReporter(output_prefix + '.csv', 12500, time=True, step=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True, volume=True, temperature=True))
# simulation.reporters.append(CheckpointReporter(output_prefix + '.chk', 12500))
simulation.reporters.append(DCDReporter(output_prefix +  '.dcd', 12500, enforcePeriodicBox=False))

# Equilibrate
print('Equilibrating...')
initial_time = time.time()
for iteration in progressbar.progressbar(range(niterations)):
    simulation.step(nsteps)
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(output_prefix + '.state.xml', 'w') as outfile:
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
