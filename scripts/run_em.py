from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

input_file = sys.argv[1]
pdb = PDBFile(input_file)
output_prefix = "/home/zhangi/choderalab/vir_collaboration/data/md_output/"
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.4) # not sure if I need to specify this (default is 7)
modeller.addSolvent(forcefield, ionicStrength=150*millimolar) # default: tip3p

system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds, hydrogenMass=4*amu)
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 4*femtoseconds) # sim temp, friction coefficient, step size

simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter(output_prefix + os.path.basename(input_file)[:-2] +  '_50ns.pdb', 12500))
simulation.reporters.append(StateDataReporter(output_prefix + os.path.basename(input_file)[:-2] + '_50ns.csv', 12500, time=True, step=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True, volume=True, temperature=True))
simulation.step(12500000)
