import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import mdtraj as md
import itertools
import numpy as np
from simtk.openmm.app import *

def identify_contacts(chain_A, chain_B, solute_reference_openmm, trajectory, frame, cutoff):
    """Identify interface residues between two chains at a given distance cutoff.
    
    Parameters
    ----------
    chain_A : int
        The index of one of the chains in the interface
    chain_B : int
        The index of the other chain in the interface
    trajectory: mdtraj.Trajectory
        The trajectory from which to choose a frame
    frame : int
        The index of the frame of the trajectory to compute distances on
    cutoff : float
        The cutoff distance definig the interface in nanometers, not inclusive

    Returns
    -------
    contacts : list of pairs (as arrays)
        The pairs of contacts, where each contact is a residue id (concatenated to the insertion code if it exists).
    """

    # Get trajectory frame
    traj = trajectory[frame]
    
    # Create dict mapping residue indexes to residue id's in PDB
    d_res = {}
    for chain in solute_reference_openmm.topology.chains():
        for res in chain.residues():
            d_res[res.index] = res.id + res.insertionCode if res.insertionCode != ' ' else res.id

    # Get residues for each chain
    residues_A = [res.index for chain in traj.topology.chains for res in chain.residues if chain.index == chain_A]
    residues_B = [res.index for chain in traj.topology.chains for res in chain.residues if chain.index == chain_B]
    
    # Enumerate all pairs of residues across chains
    pairs = list(itertools.product(residues_A, residues_B))
    
    # Calculate the distances between all residue pairs
    distance_values, distance_pairs = md.compute_contacts(traj, pairs)

    contacts = []
    for i, distance in enumerate(distance_values[0]):
        if distance < cutoff: 
            pair = distance_pairs[i]
            mapped_pair = np.array([d_res[pair[0]], d_res[pair[1]]])
            contacts.append(mapped_pair) 
    return contacts

def identify_interface_residues_4(reference_pdb, reference_withH_pdb, trajectory, cutoffs, frames, is_5udc=False):
    """Identify the interface residues for all pairs of chains involved in the interface.
    Note: This is for proteins with 4 chains.

    Parameters
    ----------
    reference_pdb : string
        The filepath for the reference PDB
    reference_withH_pdb : string
        The filepath to save the reference PDB after adding hydrogens
    trajectory : mdtraj.Trajectory
        The trajectory from which to get frames from
    cutoffs : list of floats
        The cutoff distance definig the interface in nanometers, not inclusive
    frames : list of ints
        The indexes of the frames from which to compute distances
    is_5udc : bool
        Specifies whether the pdb is a variant of 5udc, i.e. uses chains E and G for the antibody,
        instead of chains H and L.

    Returns
    -------
    rows: list of lists
        Where each list contains the pdb, cutoff, frame, and interface residues for each chain
    """

    # Load reference PDB
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    pdb = PDBFile(reference_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=7.4)
    PDBFile.writeFile(modeller.topology, modeller.positions, open(reference_withH_pdb, 'w'),  keepIds=True)
    solute_reference_mdtraj = md.load_pdb(reference_withH_pdb)
    solute_reference_openmm = PDBFile(reference_withH_pdb)

    # Load trajectory
    trajectory = md.load_dcd(trajectory, top=solute_reference_mdtraj)

    # Get interface residues
    rows = []
    for cutoff in cutoffs: 
        for frame in frames:
            if not is_5udc: # Here, H represents antibody heavy chain and L reprsents antibody light chain (not necessarily chain IDs)
                contacts_H_F = identify_contacts(0, 2, solute_reference_openmm, trajectory, frame, cutoff)
                contacts_H_X = identify_contacts(0, 3, solute_reference_openmm, trajectory, frame, cutoff)
                contacts_L_F = identify_contacts(1, 2, solute_reference_openmm, trajectory, frame, cutoff)
                contacts_L_X = identify_contacts(1, 3, solute_reference_openmm, trajectory, frame, cutoff)
            else:
                contacts_H_F = identify_contacts(4, 6, solute_reference_openmm, trajectory, frame, cutoff)
                contacts_H_X = identify_contacts(4, 7, solute_reference_openmm, trajectory, frame, cutoff)
                contacts_L_F = identify_contacts(5, 6, solute_reference_openmm, trajectory, frame, cutoff)
                contacts_L_X = identify_contacts(5, 7, solute_reference_openmm, trajectory, frame, cutoff)

            chain_H = []
            chain_L = []
            chain_F = []
            chain_X = []
                
            # Get key residues in each chain
            for pair in contacts_H_F:
                chain_H.append(pair[0])
                chain_F.append(pair[1])
            
            for pair in contacts_H_X:
                chain_H.append(pair[0])
                chain_X.append(pair[1])

            for pair in contacts_L_F:
                chain_L.append(pair[0])
                chain_F.append(pair[1])
                
            for pair in contacts_L_X:
                chain_L.append(pair[0])
                chain_X.append(pair[1])
            
            row = [os.path.basename(reference_pdb)[:-4], cutoff, frame, sorted(set(chain_H)), sorted(set(chain_L)), sorted(set(chain_F)), sorted(set(chain_X))]
            rows.append(row)
    return rows

def identify_interface_residues_12(reference_pdb, reference_withH_pdb, trajectory, cutoffs, frames):
    """Identify the interface residues for all pairs of chains involved in the interface.
    Note: This is for proteins with 12 chains.
    
    Parameters
    ----------
    reference_pdb : string
        The filepath for the reference PDB
    reference_withH_pdb : string
        The filepath to save the reference PDB after adding hydrogens
    trajectory : mdtraj.Trajectory
        The trajectory from which to get frames from
    cutoffs : list of floats
        The cutoff distance definig the interface in nanometers, not inclusive
    frames : list of ints
        The indexes of the frames from which to compute distances

    Returns
    -------
    rows: list of lists
        Where each list contains the pdb, cutoff, frame, and interface residues for each chain
    """

    # Load reference PDB
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    pdb = PDBFile(reference_pdb)
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(forcefield, pH=7.4)
    PDBFile.writeFile(modeller.topology, modeller.positions, open(reference_withH_pdb, 'w'),  keepIds=True)
    solute_reference_mdtraj = md.load_pdb(reference_withH_pdb)
    solute_reference_openmm = PDBFile(reference_withH_pdb)

    # Load trajectory
    trajectory = md.load_dcd(trajectory, top=solute_reference_mdtraj)

    # Get interface residues
    rows = []
    for cutoff in cutoffs: 
        for frame in frames:
            contacts_H_F = identify_contacts(0, 2, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_H_X = identify_contacts(0, 3, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_L_F = identify_contacts(1, 2, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_L_X = identify_contacts(1, 3, solute_reference_openmm, trajectory, frame, cutoff)
            
            contacts_B_A = identify_contacts(4, 6, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_B_Y = identify_contacts(4, 7, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_C_A = identify_contacts(5, 6, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_C_Y = identify_contacts(5, 7, solute_reference_openmm, trajectory, frame, cutoff)
            
            contacts_E_D = identify_contacts(8, 10, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_E_Z = identify_contacts(8, 11, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_G_D = identify_contacts(9, 10, solute_reference_openmm, trajectory, frame, cutoff)
            contacts_G_Z = identify_contacts(9, 11, solute_reference_openmm, trajectory, frame, cutoff)

            chain_H = []
            chain_L = []
            chain_F = []
            chain_X = []
            chain_B = []
            chain_C = []
            chain_A = []
            chain_Y = []
            chain_E = []
            chain_G = []
            chain_D = []
            chain_Z = []
                
            # Get key residues in each chain
            for pair in contacts_H_F:
                chain_H.append(pair[0])
                chain_F.append(pair[1])
            
            for pair in contacts_H_X:
                chain_H.append(pair[0])
                chain_X.append(pair[1])

            for pair in contacts_L_F:
                chain_L.append(pair[0])
                chain_F.append(pair[1])
                
            for pair in contacts_L_X:
                chain_L.append(pair[0])
                chain_X.append(pair[1])
                
            for pair in contacts_B_A:
                chain_B.append(pair[0])
                chain_A.append(pair[1])
            
            for pair in contacts_B_Y:
                chain_B.append(pair[0])
                chain_Y.append(pair[1])

            for pair in contacts_C_A:
                chain_C.append(pair[0])
                chain_A.append(pair[1])
                
            for pair in contacts_C_Y:
                chain_C.append(pair[0])
                chain_Y.append(pair[1])
            
            for pair in contacts_E_D:
                chain_E.append(pair[0])
                chain_D.append(pair[1])
            
            for pair in contacts_E_Z:
                chain_E.append(pair[0])
                chain_Z.append(pair[1])

            for pair in contacts_G_D:
                chain_G.append(pair[0])
                chain_D.append(pair[1])
                
            for pair in contacts_G_Z:
                chain_G.append(pair[0])
                chain_Z.append(pair[1])
            
            row = [os.path.basename(reference_pdb)[:-4], cutoff, frame, sorted(set(chain_H)), sorted(set(chain_L)), 
                   sorted(set(chain_F)), sorted(set(chain_X)), sorted(set(chain_B)), sorted(set(chain_C)), 
                   sorted(set(chain_A)), sorted(set(chain_Y)), sorted(set(chain_E)), sorted(set(chain_G)), 
                    sorted(set(chain_D)), sorted(set(chain_Z))]
            rows.append(row)
    return rows

# Set global parameters
combined = ['holo.nonoverlay.13', 'holo.nonoverlay.14', 'holo.overlay.4jha.9', 'holo.overlay.4jha.10', 'holo.overlay.5k6f.11', 'holo.overlay.5k6f.12']
cutoffs = [0.3, 0.5, 1.5]
frames = [0, 499, 700]
reference_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_input/renumbered/"
references = ['4jhw_final_v2_refmac1_clean_trimer_single_ab.pdb', '5udc_final_v2_refmac1_clean_single_ab.pdb', 
'4jhw_trimer_single_ab_4jha.pdb', '4jhw_trimer_single_ab_5k6f_4jha.pdb', '4jhw_trimer_single_ab_5k6f.pdb', '5udc_5k6f_single_ab.pdb']
reference_withH_prefix = reference_prefix + "addH/"
trajectory_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_output/"
trajectory_postfix = ".50ns.solute.dcd"
trajectories = ['holo.nonoverlay.13', 'holo.nonoverlay.14', 'holo.overlay.4jha.9', 
'holo.overlay.4jha.10', 'holo.overlay.5k6f.11', 'holo.overlay.5k6f.12']


# Create dataframe of interface residues for trajectories with 4 chains
rows_df_4 = []
for ref, traj in zip (references, trajectories):
    print("Trying ", os.path.basename(traj))
    reference_pdb = reference_prefix + ref
    reference_withH_pdb = reference_withH_prefix + os.path.basename(ref)[:-4] + "_withH.pdb"
    trajectory = trajectory_prefix + traj + trajectory_postfix
    if traj in combined:
        trajectory = trajectory_prefix + traj + '.50ns.combined.solute.dcd'
    rows = identify_interface_residues_4(reference_pdb, reference_withH_pdb, trajectory, cutoffs, frames)
    rows_df_4.append(rows)
rows_df_4 = [sub_row_list for row_list in rows_df_4 for sub_row_list in row_list]
pd.DataFrame(rows_df_4, columns=['trajectory_name', 'cutoff (nm)', 'frame', 'chain H', 'chain L', 'chain F', 'chain X']).to_csv("/data/chodera/zhangi/vir_collaboration/data/em_figures/df_interface_residues_4.csv", index=False)


# Create dataframe of interface residues for trajectories with 8 chains
rows_df_4 = []
for ref, traj in zip (references, trajectories):
    print("Trying ", os.path.basename(traj))
    reference_pdb = reference_prefix + ref
    reference_withH_pdb = reference_withH_prefix + os.path.basename(ref)[:-4] + "_withH.pdb"
    trajectory = trajectory_prefix + traj + trajectory_postfix
    if traj in combined:
        trajectory = trajectory_prefix + traj + '.50ns.combined.solute.dcd'
    if '5udc' in ref:
        rows = identify_interface_residues_4(reference_pdb, reference_withH_pdb, trajectory, cutoffs, frames, is_5udc=True)
    else:
        rows = identify_interface_residues_4(reference_pdb, reference_withH_pdb, trajectory, cutoffs, frames)
    rows_df_4.append(rows)
rows_df_4 = [sub_row_list for row_list in rows_df_4 for sub_row_list in row_list]
pd.DataFrame(rows_df_4, columns=['trajectory_name', 'cutoff (nm)', 'frame', 'chain Ab H', 'chain Ab L', 'chain F', 'chain X']).to_csv("/data/chodera/zhangi/vir_collaboration/data/em_figures/df_interface_residues_8.csv", index=False)

# Create dataframe of interface residues for trajectories with 12 chains (5udc full)
references = ["5udc_final_v2_refmac1_clean.pdb", "5udc_final_v2_refmac1_clean_tail.pdb", 
                "5udc_final_v2_refmac1_clean_tail_variable.pdb", "5udc_final_v2_refmac1_clean_variable.pdb", 
                "5udc_5k6f.pdb", "5udc_tail_5k6f.pdb"]
trajectories = ["holo.nonoverlay.6", "holo.nonoverlay.8", "holo.nonoverlay.10", "holo.nonoverlay.12", 
                "holo.overlay.5k6f.9", "holo.overlay.5k6f.10"]
rows_df_12 = []
for ref, traj in zip (references, trajectories):
    print("Trying ", os.path.basename(traj))
    reference_pdb = reference_prefix + ref
    reference_withH_pdb = reference_withH_prefix + os.path.basename(ref)[:-4] + "_withH.pdb"
    trajectory = trajectory_prefix + traj + trajectory_postfix
    if traj in combined:
        trajectory = trajectory_prefix + traj + '.50ns.combined.solute.dcd'
    rows = identify_interface_residues_12(reference_pdb, reference_withH_pdb, trajectory, cutoffs, frames)
    rows_df_12.append(rows)
rows_df_12 = [sub_row_list for row_list in rows_df_12 for sub_row_list in row_list]
pd.DataFrame(rows_df_12, columns=['trajectory_name', 'cutoff (nm)', 'frame', 'chain H', 'chain L', 'chain F', 'chain X', 'chain B', 'chain C', 'chain A', 'chain Y', 'chain E', 'chain G', 'chain D', 'chain Z']).to_csv("/data/chodera/zhangi/vir_collaboration/data/em_figures/df_interface_residues_12.csv", index=False)
