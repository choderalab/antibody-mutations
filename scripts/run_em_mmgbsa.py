from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import openmmtools
import pandas as pd

# Set file names
input_prefix = "/home/zhangi/choderalab/vir_collaboration/data/em_input/"
input_file_PL = input_prefix + "holo/nonoverlay/" + sys.argv[1]
input_file_P = input_prefix + "apo/" + sys.argv[2]
input_file_L = input_prefix + "apo/" + sys.argv[3]
structure_name = sys.argv[4]

# Set up forcefield
print("Loading forcefield...")
forcefield = ForceField('amber99_obc.xml', 'amber99sbildn.xml')

def addH(input_file, structure_type):
	'''Add hydrogens to structure using Modeller.

		Parameters
		----------
		input_file : string
			Path to input file of structure without Hs
		structure_type : string
			'PL', 'P', 'L'

		Returns
		-------
		modeller : simtk.openmm.app.modeller.Modeller
			Contains structure with hydrogens
	'''
	print("Adding hydrogens ", input_file)
	outfile_name = "/data/chodera/zhangi/vir_collaboration/data/mm_gbsa/" + structure_name + "_" + structure_type +  "_withH.pdb"
	pdb = PDBFile(input_file)
	modeller = Modeller(pdb.topology, pdb.positions)
	modeller.addHydrogens(forcefield, pH=7.4)
	PDBFile.writeFile(modeller.topology, modeller.positions, open(outfile_name, 'w'), keepIds=True)
	return modeller

def compute_energy(topology, positions):
	'''Run energy minimization and return initial and final energies.

		Parameters
		----------
		topology : simtk.openmm.app.topology.Topology
			Topology of structure to energy minimize
		positions: np.array
			Positions of atoms in structure to energy minimize

		Returns
		-------
		list of floats
			Initial and final energies in kcal/mol
		final_state : simtk.openmm.openmm.State
			Energy minimized final state
	'''

	# Set up system
	print('Creating OpenMM System...')
	system = forcefield.createSystem(topology, nonbondedMethod=CutoffNonPeriodic, 
		nonbondedCutoff=1.2*nanometer, constraints=HBonds, hydrogenMass=4*amu)

	# Set up integrator
	integrator = openmmtools.integrators.LangevinIntegrator(300*kelvin, 1/picosecond, 4*femtosecond)

	# Set up platform
	print("Setting up the platform...")
	platform = Platform.getPlatformByName('CUDA')
	platform.setPropertyDefaultValue('Precision', 'mixed')

	# Set up simulation 
	print("Setting up the simulation...")
	simulation = Simulation(topology, system, integrator, platform)
	simulation.context.setPositions(positions)

	# Minimize the energy
	print("Minimizing the energy...")
	initial_state = simulation.context.getState(getEnergy=True, getPositions=True)
	initial_energy = initial_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole) # before conversion: kJ/mol
	print('  initial : %8.3f kcal/mol' % (initial_energy))
	simulation.minimizeEnergy()
	final_state = simulation.context.getState(getEnergy=True, getPositions=True)
	final_energy = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
	print('  final   : %8.3f kcal/mol' % (final_energy))

	return [initial_energy, final_energy], final_state

def extract_from_PL(PL_topology, chains_to_keep, PL_positions):    
	'''Use OpenMM's Topology object to create a new topology with desired chains. 

		Parameters
		----------
		PL_topology : simtk.openmm.app.topology.Topology
			PL topology from which to extract chains
		chains_to_keep : list of ints
			Index(es) of chain(s) to keep
		PL_positions : np.array
			Positions from PL's final state

		Returns
		-------
		new_topology : simtk.openmm.app.topology.Topology
			Extracted topology of P or L 
		positions : np.array
			Sliced positions from PL_positions for P or L
	'''
	print("extracting from PL...")
	# Create new topology
	new_topology = Topology()

	# Copy residues and atoms to new topology for chains_to_keep
	d_old_to_new = {} # Key: atom in old topology, Value: atom in new topology
	for chain in PL_topology.chains():
		if chain.index in chains_to_keep:
			new_chain = new_topology.addChain(id=chain.id)
			for res in chain.residues():   
				# Copy residues and atoms
				new_res = new_topology.addResidue(res.name, new_chain, id=res.id)
				for atom in res.atoms():
					new_atom = new_topology.addAtom(atom.name, atom.element, new_res)
					d_old_to_new[atom] = new_atom

	# Make list of atoms to delete
	atoms_to_delete = []
	for res in PL_topology.residues():
		if res.chain.index not in chains_to_keep:
			for atom in res.atoms():
				atoms_to_delete.append(atom)

	# Copy bonds to new topology, except bonds involving atoms to delete
	for bond in PL_topology.bonds():
		atom_1 = bond[0]
		atom_2 = bond[1]
		if (atom_1 in atoms_to_delete) or (atom_2 in atoms_to_delete):
			continue
		atom_1_new = d_old_to_new[atom_1]
		atom_2_new = d_old_to_new[atom_2]
		new_topology.addBond(atom_1_new, atom_2_new)

	# Slice PL_positions
	atoms = [atom.index for atom in PL_topology.atoms() if atom.residue.chain.index in chains_to_keep]
	positions = [position for i, position in enumerate(PL_positions) if i in atoms]
	return new_topology, positions

# Add Hs 10x, compute energies 10x -- extract P and L from PL
all_data = []
for i in range(10):
	PL_withH = addH(input_file_PL, "PL")
	PL_data, PL_final_state = compute_energy(PL_withH.topology, PL_withH.positions)
	P_topology, P_positions = extract_from_PL(PL_withH.topology, range(2,8), PL_final_state.getPositions()) # 4jhw chain F protein
	L_topology, L_positions = extract_from_PL(PL_withH.topology, [0, 1], PL_final_state.getPositions()) # 4jhw chain F antibody
	P_data, P_final_state = compute_energy(P_topology, P_positions)
	L_data, L_final_state = compute_energy(L_topology, L_positions)
	deltaG = PL_data[1] - P_data[1] - L_data[1]
	print(deltaG)
	all_data.append([deltaG] + PL_data + P_data + L_data)
df = pd.DataFrame(all_data, columns=['deltaG', 'PL initial energy', 'PL final energy', 'P initial energy', 'P final energy', 'L initial energy', 'L final energy'])
df.to_csv("/data/chodera/zhangi/vir_collaboration/data/mm_gbsa/" + structure_name + "_initial_final_diffH_extracted.csv")
