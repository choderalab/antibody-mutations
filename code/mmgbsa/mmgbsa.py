from simtk import unit, openmm
from simtk.openmm import app
import sys
import openmmtools
import pandas as pd
import numpy as np
import subprocess

# Read in input arguments
input_file = sys.argv[1]
args = []
with open(input_file, 'r') as infile:
	args = [line.strip('\n').split("=")[1] for line in infile.readlines()]
input_file_PL_WT = args[0]
input_file_PL_mutant = args[1]
add_hydrogens = True if args[2] == 'True' else False
output_prefix = args[3]
structure_name = args[4]
protein_chains = args[5].split(",")
ligand_chains = args[6].split(",")
iterations = int(args[7]) if args[7] else 1 # Default: 1

# Check that chains were specified correctly
if not all(isinstance(chain_id, str) and len(chain_id) == 1 for chain_id in protein_chains):
	raise TypeError('Wrong format for protein chains. Should be comma-separated string of chain IDs (ex: "F,X,A,Y,D,Z")')

if not all(isinstance(chain_id, str) and len(chain_id) == 1 for chain_id in ligand_chains):
	raise TypeError('Wrong format for ligand chains. Should be comma-separated string of chain IDs (ex: "H,L")')

class MMGBSA:
	def __init__(self, input_file, structure_name):
		self._pdb = app.PDBFile(input_file)
		self._name = structure_name
		self.complexTopology = self._pdb.topology
		self.complexPositions = self._pdb.positions
		self.complexEnergy = 0
		self.proteinEnergy = 0
		self.ligandEnergy = 0
		self.minimizedState = openmm.State()
		self._pH = 7.4
		self.forcefield = app.ForceField('amber99_obc.xml', 'amber99sbildn.xml')
		self.platform = openmmtools.utils.get_fastest_platform()
		# self.platform = openmm.Platform.getPlatformByName('OpenCL')
		# self.platform.setPropertyDefaultValue('OpenCLPrecision', 'single')
		# self.platform.setPropertyDefaultValue('CudaPrecision', 'mixed')
 
	def setpH(self, pH):
		'''Set the pH at which to protonate the System, if necessary'''
		print("Setting pH...")
		if isinstance(pH, float):
			self._pH = pH
		else:
			raise TypeError("pH must be a float")

	def getpH(self):
		'''Get the pH of the System'''
		print("Getting pH...")
		return self._pH

	def addHydrogens(self, output_prefix, structure_type, iteration):
		'''Protonate the complex and update the positions and topology
		
		Parameters
		----------
		output_prefix : string
			Path to output directory in which to save protonated structure
		structure_type : string
			'WT' or 'Mutant' to be used in output file name
		iteration : int
			iteration number to use in the output file name
		'''
		print("\tAdding hydrogens...")
		modeller = app.Modeller(self.complexTopology, self.complexPositions)
		modeller.addHydrogens(self.forcefield, pH=self._pH)

		# Save protonated PDB
		print("\tWriting the protonated model to ", output_prefix + self._name + "_" + structure_type + "_" + str(iteration) + "_protonated.pdb")
		app.PDBFile.writeFile(modeller.topology, modeller.positions, open(output_prefix + self._name + "_" + structure_type + "_" + str(iteration) + "_protonated.pdb", 'w'),  keepIds=True)

		self.complexTopology = modeller.topology
		self.complexPositions = modeller.positions

	def computeEnergy(self, output_prefix, structure_type, iteration):
		'''Run energy minimization and update the final energy and final state.
		
		Parameters
		----------
		output_prefix : string
			Path to output directory in which to save minimized structure
		structure_type : string
			'WT' or 'Mutant' to be used in output file name
		iteration : int
			iteration number to use in the output file name
		Returns
		-------
		float
			Final energy in kcal/mol
		'''
		# Set up system
		system = self.forcefield.createSystem(self.complexTopology, nonbondedMethod=app.CutoffNonPeriodic, 
			nonbondedCutoff=1.2*unit.nanometer, constraints=app.HBonds, hydrogenMass=4*unit.amu, soluteDielectric=1, solventDielectric=1)

		# Set up integrator
		integrator = openmmtools.integrators.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 4*unit.femtosecond)

		# Set up simulation 
		simulation = app.Simulation(self.complexTopology, system, integrator, self.platform)
		simulation.context.setPositions(self.complexPositions)
		
		if iteration == 0:
			print("\tPlatform: ", self.platform.getName())
			print("\tPrecision: ", self.platform.getPropertyValue(simulation.context, 'Precision'))

		# Minimize the energy
		print("\tMinimizing the PL complex...")
		initial_state = simulation.context.getState(getEnergy=True, getPositions=True)
		initial_energy = initial_state.getPotentialEnergy() / unit.kilocalories_per_mole # before conversion: kJ/mol
		print('\t  initial : %8.3f kcal/mol' % (initial_energy))
		simulation.minimizeEnergy()
		final_state = simulation.context.getState(getEnergy=True, getPositions=True)
		final_energy = final_state.getPotentialEnergy() / unit.kilocalories_per_mole
		print('\t  final   : %8.3f kcal/mol' % (final_energy))
		
		# Save energy minimized PDB
		print("\tWriting the minimized model to ", output_prefix + self._name + "_" + structure_type + "_" + str(iteration) + "_minimized.pdb")
		positions_minimized = simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions()
		app.PDBFile.writeFile(self.complexTopology, positions_minimized, open(output_prefix + self._name + "_" + structure_type + "_" + str(iteration) + "_minimized.pdb", 'w'),  keepIds=True)

		self.complexEnergy = final_energy
		self.minimizedState = final_state
		return final_energy

	def extract(self, chains_to_keep, extract_protein):
		'''Extract the protein or ligand from the complex and return its energy. 

		Parameters
		----------
		chains_to_keep : list of strings
			Chain id(s) to keep
		extract_protein : boolean
			If true, extract protein. Otherwise, extract ligand.

		Returns
		-------
		float
			Final energy in kcal/mol
		'''
		
		try:
			self.minimizedState.getPositions()
		except Exception:
			print("The System has not been minimized -- you must minimize before extracting.")

		if not isinstance(extract_protein, bool):
			raise TypeError("extract_protein should be boolean indicating whether to extract the protein (True) or ligand (False)")

		## Extract topology
		if extract_protein:
			print("\tExtracting P from PL...")
		else:
			print("\tExtracting L from PL...")
		
		# Create new topology
		new_topology = app.Topology()
		
		# Copy residues and atoms to new topology for chains_to_keep
		d_old_to_new = {} # Key: atom in old topology, Value: atom in new topology
		for chain in self.complexTopology.chains():
			if chain.id in chains_to_keep:
				new_chain = new_topology.addChain(id=chain.id)
				for res in chain.residues():   
					# Copy residues and atoms
					new_res = new_topology.addResidue(res.name, new_chain, id=res.id)
					for atom in res.atoms():
						new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
						d_old_to_new[atom] = new_atom

		# Make list of atoms to delete
		atoms_to_delete = []
		for res in self.complexTopology.residues():
			if res.chain.id not in chains_to_keep:
				for atom in res.atoms():
					atoms_to_delete.append(atom)

		# Copy bonds to new topology, except bonds involving atoms to delete
		for bond in self.complexTopology.bonds():
			atom_1 = bond[0]
			atom_2 = bond[1]
			if (atom_1 in atoms_to_delete) or (atom_2 in atoms_to_delete):
				continue
			atom_1_new = d_old_to_new[atom_1]
			atom_2_new = d_old_to_new[atom_2]
			new_topology.addBond(atom_1_new, atom_2_new)

		## Extract the positions
		atoms = [atom.index for atom in self.complexTopology.atoms() if atom.residue.chain.id in chains_to_keep]
		positions = [position for i, position in enumerate(self.minimizedState.getPositions()) if i in atoms]
		
		## Extract the energy
		# Set up system
		system = self.forcefield.createSystem(new_topology, nonbondedMethod=app.CutoffNonPeriodic, 
			nonbondedCutoff=1.2*unit.nanometer, constraints=app.HBonds, hydrogenMass=4*unit.amu, soluteDielectric=1, solventDielectric=1)

		# Set up integrator
		integrator = openmmtools.integrators.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 4*unit.femtosecond)

		# Set up simulation 
		simulation = app.Simulation(new_topology, system, integrator, self.platform)
		simulation.context.setPositions(positions)

		# Get the initial energy
		initial_state = simulation.context.getState(getEnergy=True, getPositions=True)
		initial_energy = initial_state.getPotentialEnergy() / unit.kilocalories_per_mole # before conversion: kJ/mol
		print('\t  Energy : %8.3f kcal/mol' % (initial_energy))

		if extract_protein:
			self.proteinEnergy = initial_energy
		else:
			self.ligandEnergy = initial_energy

	def computeDeltaG(self):
		'''Compute MM/GBSA free energy estimate'''
		print("\tComputing MM/GBSA free energy estimate..")
		return self.complexEnergy - self.proteinEnergy - self.ligandEnergy

def main():
	all_data = []
	for i in range(iterations):
		print("\nIteration: ", i)
		print("For WT: ")
		# Calculate WT deltaG
		system = MMGBSA(input_file_PL_WT, structure_name)
		if add_hydrogens:
			system.addHydrogens(output_prefix, 'wt', i)
		system.computeEnergy(output_prefix, 'wt', i)
		system.extract(protein_chains, True)
		system.extract(ligand_chains, False)
		deltaG_wt = system.computeDeltaG()
		print('\t' + str(deltaG_wt))

		# Calculate Mutant deltaG, if mutation was specified
		if input_file_PL_mutant:
			print("For mutant: ")
			system = MMGBSA(input_file_PL_mutant, structure_name)
			if add_hydrogens:
				system.addHydrogens(output_prefix, 'mutant', i)
			system.computeEnergy(output_prefix, 'mutant', i)
			system.extract(protein_chains, True)
			system.extract(ligand_chains, False)
			deltaG_mut = system.computeDeltaG()
			print('\t' + str(deltaG_mut))
			ddG = deltaG_mut - deltaG_wt
			all_data.append([ddG, deltaG_wt, deltaG_mut])
		else:
			all_data.append([deltaG_wt])

	# Write to csv file
	if input_file_PL_mutant:
		df = pd.DataFrame(all_data, columns=['ddG', 'WT deltaG', 'Mutant deltaG'])
		df.to_csv(output_prefix + structure_name + ".csv")
	else:
		df = pd.DataFrame(all_data, columns=['WT deltaG'])
		df.to_csv(output_prefix + structure_name + ".csv")

if __name__ == "__main__":
	main()