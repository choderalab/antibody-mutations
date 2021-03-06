import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import mdtraj as md
import itertools
import numpy as np
from simtk.openmm.app import *

def plot_interface_rmsd_4(df, reference_withH_pdb, trajectory, frame, outfile_prefix):
    """Align on each of the components and generate component RMSD plots for each of the alignments.
    Note: this is for trajectories with four chains.
    
    Parameters
    ----------
    df : pd.DataFrame
        Contains the interface residues for each chain in each trajectory
    reference_withH_pdb : string
        The filepath to the reference pdb with hydrogens added
    trajectory: string
        The filepath to the trajectory
    frame : int
        The frame number in trajectory on which to compute distances upon
    outfile_prefix : string
        The output file prefix containing directories to save the rmsd plots
    """

    # Load reference PDB
    solute_reference = md.load_pdb(reference_withH_pdb)

    # Load trajectory
    trajectory = md.load_dcd(trajectory, top=solute_reference)
    traj_name = os.path.basename(reference_withH_pdb[:-10])
    
    # Image molecules against first chain
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[0:1])
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[1:2])

    # Get interface residues for cutoff = 5 angstroms
    row = df[(df['trajectory_name'] == traj_name) & (df['cutoff (nm)'] == 0.5) & (df['frame'] == frame)]
    chain_H = row['chain H']
    chain_L = row['chain L']
    chain_F = row['chain F']
    chain_X = row['chain X']

    # Create string for interface selection query and select interface atoms
    chain_H_str = ' or resSeq '.join(chain_H.values[0][1:-1].replace("'", "").split(", "))
    chain_L_str = ' or resSeq '.join(chain_L.values[0][1:-1].replace("'", "").split(", "))
    chain_F_str = ' or resSeq '.join(chain_F.values[0][1:-1].replace("'", "").split(", "))
    chain_X_str = ' or resSeq '.join(chain_X.values[0][1:-1].replace("'", "").split(", "))
    selection_str = '(chainid 0 and (resSeq ' + chain_H_str + '))' + ' or (chainid 1 and (resSeq ' + chain_L_str + '))' + ' or (chainid 2 and (resSeq ' + chain_F_str + '))' + ' or (chainid 3 and (resSeq ' + chain_X_str + '))'
    interface_atoms_b_5 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    interface_atoms_h_5 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Get interface residues for cutoff = 15 angstroms
    row = df[(df['trajectory_name'] == traj_name) & (df['cutoff (nm)'] == 1.5) & (df['frame'] == frame)]
    chain_H = row['chain H']
    chain_L = row['chain L']
    chain_F = row['chain F']
    chain_X = row['chain X']

    # Create string for interface selection query and select interface atoms
    chain_H_str = ' or resSeq '.join(chain_H.values[0][1:-1].replace("'", "").split(", "))
    chain_L_str = ' or resSeq '.join(chain_L.values[0][1:-1].replace("'", "").split(", "))
    chain_F_str = ' or resSeq '.join(chain_F.values[0][1:-1].replace("'", "").split(", "))
    chain_X_str = ' or resSeq '.join(chain_X.values[0][1:-1].replace("'", "").split(", "))
    selection_str = '(chainid 0 and (resSeq ' + chain_H_str + '))' + ' or (chainid 1 and (resSeq ' + chain_L_str + '))' + ' or (chainid 2 and (resSeq ' + chain_F_str + '))' + ' or (chainid 3 and (resSeq ' + chain_X_str + '))'
    interface_atoms_b_15 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    interface_atoms_h_15 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Create string for F protein selection query and select interface F protein atoms
    selection_str = '(chainid 2) or (chainid 3)' 
    f_atoms_b = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    f_atoms_h = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Create string for selection query
    selection_str = '(chainid 0) or (chainid 1)' 
    antibody_atoms_b = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    antibody_atoms_h = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Align on interface (5 angstroms), compute and plot RMSDs
    trajectory_i_5 = trajectory.superpose(solute_reference, frame=frame, atom_indices=interface_atoms_b_5, ref_atom_indices=interface_atoms_b_5)
    compute_and_plot(trajectory_i_5, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'interface (5)')
    
    # Align on interface (15 angstroms), compute and plot RMSDs
    trajectory_i_15 = trajectory.superpose(solute_reference, frame=frame, atom_indices=interface_atoms_b_15, ref_atom_indices=interface_atoms_b_15)
    compute_and_plot(trajectory_i_15, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'interface (15)')

    # Align to F protein, compute and plot RMSDs
    trajectory_f = trajectory.superpose(solute_reference, frame=frame, atom_indices=f_atoms_b, ref_atom_indices=f_atoms_b)
    compute_and_plot(trajectory_f, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'F protein')
    
    # Align to antibody, compute and plot RMSDs
    trajectory_a = trajectory.superpose(solute_reference, frame=frame, atom_indices=antibody_atoms_b, ref_atom_indices=antibody_atoms_b)
    compute_and_plot(trajectory_a, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'antibody')

def plot_interface_rmsd_8(df, reference_withH_pdb, trajectory, frame, outfile_prefix, is_5udc=False):
    """Align on each of the components and generate component RMSD plots for each of the alignments.
    Note: this is for trajectories with eight chains.
    
    Parameters
    ----------
    df : pd.DataFrame
        Contains the interface residues for each chain in each trajectory
    reference_withH_pdb : string
        The filepath to the reference pdb with hydrogens added
    trajectory: string
        The filepath to the trajectory
    frame : int
        The frame number in trajectory on which to compute distances upon
    outfile_prefix : string
        The output file prefix containing directories to save the rmsd plots
    is_5udc : bool
        Specifies whether the pdb is a variant of 5udc, i.e. uses chains E and G for the antibody,
        instead of chains H and L.
    """
    # Load reference PDB
    solute_reference = md.load_pdb(reference_withH_pdb)

    # Load trajectory
    trajectory = md.load_dcd(trajectory, top=solute_reference)
    traj_name = os.path.basename(reference_withH_pdb[:-10])
    
    # Image molecules against first chain
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[0:1])
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[1:2])
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[2:3])

    # Get interface residues for cutoff = 5 angstroms
    row = df[(df['trajectory_name'] == traj_name) & (df['cutoff (nm)'] == 0.5) & (df['frame'] == frame)]
    
    # Antibody chains
    chain_H = row['chain Ab H']
    chain_L = row['chain Ab L']

    # F protein chains
    chain_F = row['chain F']
    chain_X = row['chain X']

    # Create string for interface selection query and select interface atoms
    chain_H_str = ' or resSeq '.join(chain_H.values[0][1:-1].replace("'", "").split(", "))
    chain_L_str = ' or resSeq '.join(chain_L.values[0][1:-1].replace("'", "").split(", "))
    chain_F_str = ' or resSeq '.join(chain_F.values[0][1:-1].replace("'", "").split(", "))
    chain_X_str = ' or resSeq '.join(chain_X.values[0][1:-1].replace("'", "").split(", "))

    if not is_5udc:
        selection_str = '(chainid 0 and (resSeq ' + chain_H_str + '))' + ' or (chainid 1 and (resSeq ' + chain_L_str + '))' \
        + ' or (chainid 2 and (resSeq ' + chain_F_str + '))' + ' or (chainid 3 and (resSeq ' + chain_X_str + '))'
    else:
        selection_str = '(chainid 4 and (resSeq ' + chain_H_str + '))' + ' or (chainid 5 and (resSeq ' + chain_L_str + '))' \
        + ' or (chainid 0 and (resSeq ' + chain_F_str + '))' + ' or (chainid 1 and (resSeq ' + chain_X_str + '))'

    interface_atoms_b_5 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    interface_atoms_h_5 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Get interface residues for cutoff = 15 angstroms
    row = df[(df['trajectory_name'] == traj_name) & (df['cutoff (nm)'] == 1.5) & (df['frame'] == frame)]
    
    # Antibody chains
    chain_H = row['chain Ab H']
    chain_L = row['chain Ab L']

    # F protein chains
    chain_F = row['chain F']
    chain_X = row['chain X']

    # Create string for interface selection query and select interface atoms
    chain_H_str = ' or resSeq '.join(chain_H.values[0][1:-1].replace("'", "").split(", "))
    chain_L_str = ' or resSeq '.join(chain_L.values[0][1:-1].replace("'", "").split(", "))
    chain_F_str = ' or resSeq '.join(chain_F.values[0][1:-1].replace("'", "").split(", "))
    chain_X_str = ' or resSeq '.join(chain_X.values[0][1:-1].replace("'", "").split(", "))
    
    if not is_5udc:
        selection_str = '(chainid 0 and (resSeq ' + chain_H_str + '))' + ' or (chainid 1 and (resSeq ' + chain_L_str + '))' \
        + ' or (chainid 2 and (resSeq ' + chain_F_str + '))' + ' or (chainid 3 and (resSeq ' + chain_X_str + '))'
    else:
        selection_str = '(chainid 4 and (resSeq ' + chain_H_str + '))' + ' or (chainid 5 and (resSeq ' + chain_L_str + '))' \
        + ' or (chainid 0 and (resSeq ' + chain_F_str + '))' + ' or (chainid 1 and (resSeq ' + chain_X_str + '))'

    interface_atoms_b_15 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    interface_atoms_h_15 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Create string for F protein selection query and select interface F protein atoms
    if not is_5udc:
        selection_str = '(chainid 2) or (chainid 3) or (chainid 4) or (chainid 5) or (chainid 6) or (chainid 7)' 
    else:
        selection_str = '(chainid 0) or (chainid 1) or (chainid 2) or (chainid 3) or (chainid 6) or (chainid 7)' 

    f_atoms_b = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    f_atoms_h = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Create string for selection query antibody selection query and select interface antibody atoms
    if not is_5udc:
        selection_str = '(chainid 0) or (chainid 1)' 
    else:
        selection_str = '(chainid 4) or (chainid 5)' 
    antibody_atoms_b = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    antibody_atoms_h = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Align on interface (5 angstroms), compute and plot RMSDs
    trajectory_i_5 = trajectory.superpose(solute_reference, frame=frame, atom_indices=interface_atoms_b_5, ref_atom_indices=interface_atoms_b_5)
    compute_and_plot(trajectory_i_5, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'interface (5)')
    
    # Align on interface (15 angstroms), compute and plot RMSDs
    trajectory_i_15 = trajectory.superpose(solute_reference, frame=frame, atom_indices=interface_atoms_b_15, ref_atom_indices=interface_atoms_b_15)
    compute_and_plot(trajectory_i_15, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'interface (15)')

    # Align to F protein, compute and plot RMSDs
    trajectory_f = trajectory.superpose(solute_reference, frame=frame, atom_indices=f_atoms_b, ref_atom_indices=f_atoms_b)
    compute_and_plot(trajectory_f, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'F protein')
    
    # Align to antibody, compute and plot RMSDs
    trajectory_a = trajectory.superpose(solute_reference, frame=frame, atom_indices=antibody_atoms_b, ref_atom_indices=antibody_atoms_b)
    compute_and_plot(trajectory_a, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'antibody')


def plot_interface_rmsd_12(df, reference_withH_pdb, trajectory, frame, outfile_prefix):
    """Align on each of the components and generate component RMSD plots for each of the alignments.
    Note: this is for trajectories with twelve chains.
    
    Parameters
    ----------
    df : pd.DataFrame
        Contains the interface residues for each chain in each trajectory
    reference_withH_pdb : string
        The filepath to the reference pdb with hydrogens added
    trajectory: string
        The filepath to the trajectory
    frame : int
        The frame number in trajectory on which to compute distances upon
    outfile_prefix : string
        The output file prefix containing directories to save the rmsd plots
    """
    # Load reference PDB
    solute_reference = md.load_pdb(reference_withH_pdb)

    # Load trajectory
    trajectory = md.load_dcd(trajectory, top=solute_reference)
    traj_name = os.path.basename(reference_withH_pdb[:-10])
    
    # Image molecules against first chain
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[0:1])
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[1:2])

    # Get interface residues for cutoff = 5 angstroms
    row = df[(df['trajectory_name'] == traj_name) & (df['cutoff (nm)'] == 0.5) & (df['frame'] == frame)]
    chain_H = row['chain H']
    chain_L = row['chain L']
    chain_F = row['chain F']
    chain_X = row['chain X']
    chain_B = row['chain B']
    chain_C = row['chain C']
    chain_A = row['chain A']
    chain_Y = row['chain Y']
    chain_E = row['chain E']
    chain_G = row['chain G']
    chain_D = row['chain D']
    chain_Z = row['chain Z']

    # Create string for interface selection query and select interface atoms
    chain_H_str = ' or resSeq '.join(chain_H.values[0][1:-1].replace("'", "").split(", "))
    chain_L_str = ' or resSeq '.join(chain_L.values[0][1:-1].replace("'", "").split(", "))
    chain_F_str = ' or resSeq '.join(chain_F.values[0][1:-1].replace("'", "").split(", "))
    chain_X_str = ' or resSeq '.join(chain_X.values[0][1:-1].replace("'", "").split(", "))
    chain_B_str = ' or resSeq '.join(chain_B.values[0][1:-1].replace("'", "").split(", "))
    chain_C_str = ' or resSeq '.join(chain_C.values[0][1:-1].replace("'", "").split(", "))
    chain_A_str = ' or resSeq '.join(chain_A.values[0][1:-1].replace("'", "").split(", "))
    chain_Y_str = ' or resSeq '.join(chain_Y.values[0][1:-1].replace("'", "").split(", "))
    chain_E_str = ' or resSeq '.join(chain_E.values[0][1:-1].replace("'", "").split(", "))
    chain_G_str = ' or resSeq '.join(chain_G.values[0][1:-1].replace("'", "").split(", "))
    chain_D_str = ' or resSeq '.join(chain_D.values[0][1:-1].replace("'", "").split(", "))
    chain_Z_str = ' or resSeq '.join(chain_Z.values[0][1:-1].replace("'", "").split(", "))
    selection_str = '(chainid 0 and (resSeq ' + chain_H_str + '))' + ' or (chainid 1 and (resSeq ' + chain_L_str + '))' \
    + ' or (chainid 2 and (resSeq ' + chain_F_str + '))' + ' or (chainid 3 and (resSeq ' + chain_X_str + '))' \
    + ' or (chainid 4 and (resSeq ' + chain_B_str + '))' + ' or (chainid 5 and (resSeq ' + chain_C_str + '))' \
    + ' or (chainid 6 and (resSeq ' + chain_A_str + '))' + ' or (chainid 7 and (resSeq ' + chain_Y_str + '))' \
    + ' or (chainid 8 and (resSeq ' + chain_E_str + '))' + ' or (chainid 9 and (resSeq ' + chain_G_str + '))' \
    + ' or (chainid 10 and (resSeq ' + chain_D_str + '))' + ' or (chainid 11 and (resSeq ' + chain_Z_str + '))'
    interface_atoms_b_5 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    interface_atoms_h_5 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

   # Get interface residues for cutoff = 15 angstroms
    row = df[(df['trajectory_name'] == traj_name) & (df['cutoff (nm)'] == 1.5) & (df['frame'] == frame)]
    chain_H = row['chain H']
    chain_L = row['chain L']
    chain_F = row['chain F']
    chain_X = row['chain X']
    chain_B = row['chain B']
    chain_C = row['chain C']
    chain_A = row['chain A']
    chain_Y = row['chain Y']
    chain_E = row['chain E']
    chain_G = row['chain G']
    chain_D = row['chain D']
    chain_Z = row['chain Z']

    # Create string for interface selection query and select interface atoms
    chain_H_str = ' or resSeq '.join(chain_H.values[0][1:-1].replace("'", "").split(", "))
    chain_L_str = ' or resSeq '.join(chain_L.values[0][1:-1].replace("'", "").split(", "))
    chain_F_str = ' or resSeq '.join(chain_F.values[0][1:-1].replace("'", "").split(", "))
    chain_X_str = ' or resSeq '.join(chain_X.values[0][1:-1].replace("'", "").split(", "))
    chain_B_str = ' or resSeq '.join(chain_B.values[0][1:-1].replace("'", "").split(", "))
    chain_C_str = ' or resSeq '.join(chain_C.values[0][1:-1].replace("'", "").split(", "))
    chain_A_str = ' or resSeq '.join(chain_A.values[0][1:-1].replace("'", "").split(", "))
    chain_Y_str = ' or resSeq '.join(chain_Y.values[0][1:-1].replace("'", "").split(", "))
    chain_E_str = ' or resSeq '.join(chain_E.values[0][1:-1].replace("'", "").split(", "))
    chain_G_str = ' or resSeq '.join(chain_G.values[0][1:-1].replace("'", "").split(", "))
    chain_D_str = ' or resSeq '.join(chain_D.values[0][1:-1].replace("'", "").split(", "))
    chain_Z_str = ' or resSeq '.join(chain_Z.values[0][1:-1].replace("'", "").split(", "))
    selection_str = '(chainid 0 and (resSeq ' + chain_H_str + '))' + ' or (chainid 1 and (resSeq ' + chain_L_str + '))' \
    + ' or (chainid 2 and (resSeq ' + chain_F_str + '))' + ' or (chainid 3 and (resSeq ' + chain_X_str + '))' \
    + ' or (chainid 4 and (resSeq ' + chain_B_str + '))' + ' or (chainid 5 and (resSeq ' + chain_C_str + '))' \
    + ' or (chainid 6 and (resSeq ' + chain_A_str + '))' + ' or (chainid 7 and (resSeq ' + chain_Y_str + '))' \
    + ' or (chainid 8 and (resSeq ' + chain_E_str + '))' + ' or (chainid 9 and (resSeq ' + chain_G_str + '))' \
    + ' or (chainid 10 and (resSeq ' + chain_D_str + '))' + ' or (chainid 11 and (resSeq ' + chain_Z_str + '))'
    interface_atoms_b_15 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    interface_atoms_h_15 = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Create string for F protein selection query and select interface F protein atoms
    selection_str = '(chainid 2) or (chainid 3) or (chainid 6) or (chainid 7) or (chainid 10) or (chainid 11)' 
    f_atoms_b = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    f_atoms_h = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Create string for selection query
    selection_str = '(chainid 0) or (chainid 1) or (chainid 4) or (chainid 5) or (chainid 8) or (chainid 9)' 
    antibody_atoms_b = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select('backbone')
    antibody_atoms_h = trajectory.atom_slice(trajectory.topology.select(selection_str)).topology.select_atom_indices('heavy')

    # Align on interface (5 angstroms), compute and plot RMSDs
    trajectory_i_5 = trajectory.superpose(solute_reference, frame=frame, atom_indices=interface_atoms_b_5, ref_atom_indices=interface_atoms_b_5)
    compute_and_plot(trajectory_i_5, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
    	interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
    	antibody_atoms_h, 'interface (5)')
    
    # Align on interface (15 angstroms), compute and plot RMSDs
    trajectory_i_15 = trajectory.superpose(solute_reference, frame=frame, atom_indices=interface_atoms_b_15, ref_atom_indices=interface_atoms_b_15)
    compute_and_plot(trajectory_i_15, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'interface (15)')

    # Align to F protein, compute and plot RMSDs
    trajectory_f = trajectory.superpose(solute_reference, frame=frame, atom_indices=f_atoms_b, ref_atom_indices=f_atoms_b)
    compute_and_plot(trajectory_f, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'F protein')
    
    # Align to antibody, compute and plot RMSDs
    trajectory_a = trajectory.superpose(solute_reference, frame=frame, atom_indices=antibody_atoms_b, ref_atom_indices=antibody_atoms_b)
    compute_and_plot(trajectory_a, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
        interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
        antibody_atoms_h, 'antibody')

def compute_and_plot(aligned_trajectory, traj_name, outfile_prefix, solute_reference, interface_atoms_b_5, 
    interface_atoms_h_5, interface_atoms_b_15, interface_atoms_h_15, f_atoms_b, f_atoms_h, antibody_atoms_b, 
    antibody_atoms_h, alignment_type):
    """Computes and plot interface RMSD given interface atoms at different cutoffs, aligned on different components.
    
    Parameters
    ----------
    aligned_trajectory : mdtraj.Trajectory
        The trajectory aligned on certain components
    traj_name : string
        The name of the trajetory
    outfile_prefix : string
        The filepath at which to store the plots
    solute_reference : mdtraj.Trajectory
        The solute PDB used as reference for calculating RMSD
    interface_atoms_b_5 : list of ints
        The indices of backbone atoms in the interface as defined by a 5 angstrom cutoff
    interface_atoms_h_5 : list of ints
        The indices of heavy atoms in the interface as defined by a 5 angstrom cutoff
    interface_atoms_b_15 : list of ints
        The indices of backbone atoms in the interface as defined by a 15 angstrom cutoff
    interface_atoms_h_15 : list of ints
        The indices of heavy atoms in the interface as defined by a 15 angstrom cutoff
    f_atoms_b : list of ints
        The indices of backbone atoms in the F protein
    f_atoms_h : list of ints
        The indices of heavy atoms in the F protein
    antibody_atoms_b : list of ints
        The indices of backbone atoms in the antibody
    antibody_atoms_h : list of ints
        The indices of heavy atoms in the F antibody
    alignment_type : string
        Indicates the type of alignment (options: 'antibody', 'F protein', 'interface (15)', 'interface (5)')
    """

    # Compute interface RMSD manually without additional alignment
    rmsds_i_b_5 = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, interface_atoms_b_5, :] - solute_reference.xyz[:, interface_atoms_b_5, :])**2, axis=(1,2)))
    rmsds_i_h_5 = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, interface_atoms_h_5, :] - solute_reference.xyz[:, interface_atoms_h_5, :])**2, axis=(1,2)))

    # Compute interface RMSD manually without additional alignment
    rmsds_i_b_15 = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, interface_atoms_b_15, :] - solute_reference.xyz[:, interface_atoms_b_15, :])**2, axis=(1,2)))
    rmsds_i_h_15 = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, interface_atoms_h_15, :] - solute_reference.xyz[:, interface_atoms_h_15, :])**2, axis=(1,2)))

    # Compute the F protein RMSD manually without additional alignment
    rmsds_f_b = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, f_atoms_b, :] - solute_reference.xyz[:, f_atoms_b, :])**2, axis=(1,2)))
    rmsds_f_h = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, f_atoms_h, :] - solute_reference.xyz[:, f_atoms_h, :])**2, axis=(1,2)))

    # Compute the antibody RMSD manually without additional alignment
    rmsds_antibody_b = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, antibody_atoms_b, :] - solute_reference.xyz[:, antibody_atoms_b, :])**2, axis=(1,2)))
    rmsds_antibody_h = np.sqrt(3*np.mean((aligned_trajectory.xyz[:, antibody_atoms_h, :] - solute_reference.xyz[:, antibody_atoms_h, :])**2, axis=(1,2)))

    # Plot RMSD
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_i_b_5 * 10, color='dodgerblue', label='interface (5) backbone atoms', ax=ax)
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_i_h_5 * 10, color='cyan', label='interface (5) heavy atoms', ax=ax)
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_i_b_15 * 10, color='forestgreen', label='interface (15) backbone atoms', ax=ax)
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_i_h_15 * 10, color='lime', label='interface (15) heavy atoms', ax=ax)
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_f_b * 10, color='r', label='f backbone atoms', ax=ax)
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_f_h * 10, color='lightpink', label='f heavy atoms', ax=ax)
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_antibody_b * 10, color='darkorange', label='antibody backbone atoms', ax=ax)
    sns.lineplot(aligned_trajectory.time * 50 / 1000, rmsds_antibody_h * 10, color='gold', label='antibody heavy atoms', ax=ax)
    plt.legend()
    ax.set(xlabel='time (ns)', ylabel='RMSD (angstroms)', title=traj_name + "\naligned on " + alignment_type)
    plt.savefig(outfile_prefix + traj_name + "_" + alignment_type + "_rmsd.png", dpi=500) 

# Set global parameters
frame = 0
outfile_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_figures/"
df_4_orig = pd.read_csv("/data/chodera/zhangi/vir_collaboration/data/em_figures/df_interface_residues_4.csv")
df_8 = pd.read_csv("/data/chodera/zhangi/vir_collaboration/data/em_figures/df_interface_residues_8.csv")
df_12 = pd.read_csv("/data/chodera/zhangi/vir_collaboration/data/em_figures/df_interface_residues_12_combined.csv")


combined = ['holo.nonoverlay.5', 'holo.nonoverlay.6', 'holo.overlay.4jha.1', 'holo.overlay.4jha.5', 
'holo.overlay.5k6f.1', 'holo.overlay.5k6f.5', 'holo.nonoverlay.13', 'holo.nonoverlay.14', 'holo.overlay.4jha.9', 
'holo.overlay.4jha.10', 'holo.overlay.5k6f.9', 'holo.overlay.5k6f.10' 'holo.overlay.5k6f.11', 'holo.overlay.5k6f.12'] # Trajectories that needed to be restarted and combined
references = ["4jhw_final_v2_refmac1_clean", "4jhw_final_v2_refmac1_clean_tail", 
                "4jhw_final_v2_refmac1_clean_tail_variable", "4jhw_final_v2_refmac1_clean_variable",
              "5udc_final_v2_refmac1_clean_monomer", "5udc_final_v2_refmac1_clean_tail_monomer", 
              "5udc_final_v2_refmac1_clean_tail_variable_monomer", "5udc_final_v2_refmac1_clean_variable_monomer",
               "4jhw_5k6f_4jha", "4jhw_5k6f_tail_4jha", "4jhw_variable_5k6f_4jha", "4jhw_variable_5k6f_tail_4jha", 
              "4jhw_5k6f", "4jhw_5k6f_tail", "4jhw_variable_5k6f", "4jhw_variable_5k6f_tail", 
              "5udc_monomer_5k6f", "5udc_tail_monomer_5k6f",
              "5udc_tail_variable_monomer_5k6f", "5udc_variable_monomer_5k6f", "4jhw_final_v2_refmac1_clean_trimer_single_ab.pdb", 
              "5udc_final_v2_refmac1_clean_single_ab.pdb", "4jhw_trimer_single_ab_4jha.pdb", "4jhw_trimer_single_ab_5k6f_4jha.pdb", 
              "4jhw_trimer_single_ab_5k6f.pdb", "5udc_5k6f_single_ab.pdb"]
reference_withH_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_input/renumbered/addH/"
reference_withH_postfix = "_withH.pdb"
trajectory_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_output/"
trajectory_postfix = ".50ns.solute.dcd"
trajectories = ["holo.nonoverlay.1", "holo.nonoverlay.2", "holo.nonoverlay.3", "holo.nonoverlay.4", 
                "holo.nonoverlay.5", "holo.nonoverlay.7", "holo.nonoverlay.9", "holo.nonoverlay.11",
               "holo.overlay.4jha.1",  "holo.overlay.4jha.2", "holo.overlay.4jha.3", "holo.overlay.4jha.4", 
               "holo.overlay.5k6f.1", "holo.overlay.5k6f.2", "holo.overlay.5k6f.3", "holo.overlay.5k6f.4", 
               "holo.overlay.5k6f.5", "holo.overlay.5k6f.6", "holo.overlay.5k6f.7", "holo.overlay.5k6f.8", 
               "holo.nonoverlay.13", "holo.nonoverlay.14", "holo.overlay.4jha.9", "holo.overlay.4jha.10", 
               "holo.overlay.5k6f.11", "holo.overlay.5k6f.12"]

# Plot interface/component RMSDs for 4 chain trajectories
for ref, traj in zip(references, trajectories):
    print("Trying: ", traj)
    reference_withH_pdb = reference_withH_prefix + ref + reference_withH_postfix
    trajectory = trajectory_prefix + traj + trajectory_postfix
    if traj in combined:
        trajectory = trajectory_prefix + traj + '.50ns.combined.solute.dcd'
    plot_interface_rmsd_4(df_4, reference_withH_pdb, trajectory, frame, outfile_prefix)

# Plot interface/component RMSDs for 8 chain trajectories
references = ["4jhw_final_v2_refmac1_clean_trimer_single_ab", "5udc_final_v2_refmac1_clean_single_ab", 
                "4jhw_trimer_single_ab_4jha", "4jhw_trimer_single_ab_5k6f_4jha", 
                "4jhw_trimer_single_ab_5k6f", "5udc_5k6f_single_ab_5k6f"]
trajectories = ["holo.nonoverlay.13", "holo.nonoverlay.14", "holo.overlay.4jha.9", "holo.overlay.4jha.10",
"holo.overlay.5k6f.11", "holo.overlay.5k6f.12"]

for ref, traj in zip(references, trajectories):
    print("Trying: ", traj)
    reference_withH_pdb = reference_withH_prefix + ref + reference_withH_postfix
    trajectory = trajectory_prefix + traj + trajectory_postfix
    if traj in combined:
        trajectory = trajectory_prefix + traj + '.50ns.combined.solute.dcd'
    if '5udc' in ref:
        plot_interface_rmsd_8(df_8, reference_withH_pdb, trajectory, frame, outfile_prefix, is_5udc=True)
    else:
        plot_interface_rmsd_8(df_8, reference_withH_pdb, trajectory, frame, outfile_prefix)

# Plot interface/component RMSDs for 12 chain trajectories
references = ["5udc_final_v2_refmac1_clean", "5udc_final_v2_refmac1_clean_tail", 
                "5udc_final_v2_refmac1_clean_tail_variable", "5udc_final_v2_refmac1_clean_variable", 
                "5udc_5k6f", "5udc_tail_5k6f"]
trajectories = ["holo.nonoverlay.6", "holo.nonoverlay.8", "holo.nonoverlay.10", "holo.nonoverlay.12",
"holo.overlay.5k6f.9", "holo.overlay.5k6f.10"]

for ref, traj in zip(references, trajectories):
    print("Trying: ", traj)
    reference_withH_pdb = reference_withH_prefix + ref + reference_withH_postfix
    trajectory = trajectory_prefix + traj + trajectory_postfix
    if traj in combined:
        trajectory = trajectory_prefix + traj + '.50ns.combined.solute.dcd'
    plot_interface_rmsd_12(df_12, reference_withH_pdb, trajectory, frame, outfile_prefix)
