import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import mdtraj as md

# Description: Given the state data file for the simulation, plot all variables against time and save plots as one figure.
# state_data_file : (str) file path to state data from simulation
# title : (str) desired title for figure
def plot_all(state_data_file, title):
    print("Plotting state data...")
    data = pd.read_csv(state_data_file)
    d_labels = {"time": "Time (ps)", 
                "PE":'Potential Energy (kJ/mole)', 
                "KE": "Kinetic Energy (kJ/mole)", 
                "T": "Temperature (K)", 
                "V": "Box Volume (nm^3)"}
    fig, ax = plt.subplots(2,2, figsize=(15,10))
    sns.lineplot(x=d_labels['time'], y=d_labels['PE'], data=data, ax=ax[0][0])
    sns.lineplot(x=d_labels['time'], y=d_labels['KE'], data=data, ax=ax[1][0])
    sns.lineplot(x=d_labels['time'], y=d_labels['T'], data=data,  ax=ax[0][1])
    sns.scatterplot(x=d_labels['time'], y=d_labels['V'], data=data,  ax=ax[1][1])
    fig.suptitle(title, position=(0.5, 0.93), fontsize=20)
    plt.savefig(outfile_prefix + title + "_state_data.png", dpi=500)

# Description: Given the trajectory file and minimized system pdb for the simulation, plot RMSD against each frame and save figure.
# trajectory_file : (str) file path to trajectory from simulation
# reference_pdb : (str) file path to minimized system pdb from simulation
# title : (str) desired title for figure
def plot_rmsd(trajectory_file, reference_pdb, title):
    # Load reference pdb and extract solvent
    print("Loading reference PDB...")
    reference = md.load_pdb(reference_pdb)
    solute_reference = reference.remove_solvent()
    
    # Load trajectory, center, and superpose it to the reference
    print("Loading trajectory...")
    trajectory = md.load_dcd(trajectory_file, top=solute_reference)

    # Image molecules against first chain
    print("Imaging molecules...")
    trajectory = trajectory.image_molecules(anchor_molecules=solute_reference.topology.find_molecules()[0:1])

    # Calculated RMSD for backbone atoms vs. heavy atoms
    print("Calculating RMSD...")
    backbone_atoms = trajectory.topology.select('backbone')
    heavy_atoms = trajectory.topology.select_atom_indices('heavy')
    rmsds = md.rmsd(trajectory, solute_reference, 0, atom_indices=backbone_atoms)
    rmsds_h = md.rmsd(trajectory, solute_reference, 0, atom_indices=heavy_atoms)
    
    # Plot
    print("Plotting RMSD...")
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    sns.lineplot(trajectory.time * 50, rmsds * 10, color='r', label='backbone atoms', ax=ax)
    sns.lineplot(trajectory.time * 50, rmsds_h * 10, color='b', label='heavy atoms', ax=ax)
    plt.legend()
    ax.set(xlabel='time (ps)', ylabel='RMSD (angstroms)', title=title)
    plt.show()
    plt.savefig(outfile_prefix + title + "_rmsd.png", dpi=500)

# Create dict where key : trajectory file name, value : desired title for the plots
d = {'apo.1': '4jha_clean', 'apo.2': '4jha_clean_variable', 'apo.3': '5k6f_splitchain_capped', 
    'apo.4': '5k6f_splitchain_capped_tail', 'apo.5': '4jhw_antibody', 
    'apo.6': '4jhw_antibody_variable', 'apo.7': '5udc_antibody',
    'apo.8': '5udc_antibody_variable','holo.nonoverlay.1': '4jhw_rerefined_clean', 
    'holo.nonoverlay.2': '4jhw_rerefined_clean_tail', 'holo.nonoverlay.3': '4jhw_rerefined_clean_tail_variable', 
    'holo.nonoverlay.4': '4jhw_rerefined_clean_variable', 'holo.nonoverlay.5': '5udc_rerefined_clean_monomer',
    'holo.nonoverlay.6':'5udc_rerefined_clean', 'holo.nonoverlay.7': '5udc_rerefined_clean_tail_monomer', 
    'holo.nonoverlay.8': '5udc_rerefined_clean_tail', 'holo.nonoverlay.9':'5udc_rerefined_clean_tail_variable_monomer', 
    'holo.nonoverlay.10':'5udc_rerefined_clean_tail_variable', 'holo.nonoverlay.11':'5udc_clean_variable_monomer', 
    'holo.nonoverlay.12': '5udc_rerefined_clean_variable', 'holo.overlay.4jha.1':'4jhw_5k6f_4jha',
    'holo.overlay.4jha.2': '4jhw_5k6f_tail_4jha', 'holo.overlay.4jha.3':'4jhw_variable_5k6f_4jha', 
    'holo.overlay.4jha.4':'4jhw_variable_5k6f_tail_4jha', 'holo.overlay.4jha.5': '5udc_monomer_5k6f_4jha', 
    'holo.overlay.4jha.6': '5udc_tail_monomer_5k6f_4jha', 'holo.overlay.4jha.7': '5udc_tail_variable_monomer_5k6f_4jha',
    'holo.overlay.4jha.8': '5udc_variable_monomer_5k6f_4jha', 'holo.overlay.5k6f.1': '4jhw_5k6f', 
    'holo.overlay.5k6f.2': '4jhw_5k6f_tail', 'holo.overlay.5k6f.3': '4jhw_variable_5k6f', 
    'holo.overlay.5k6f.4': '4jhw_variable_5k6f_tail', 'holo.overlay.5k6f.5': '5udc_monomer_5k6f', 
    'holo.overlay.5k6f.6': '5udc_tail_monomer_5k6f', 'holo.overlay.5k6f.7': '5udc_tail_variable_monomer_5k6f',
    'holo.overlay.5k6f.8': '5udc_variable_monomer_5k6f', 'apo.9': '5udc_antibody', 
    'apo.10': '5udc_variable_antibody', 'apo.11': '4jhw_f', 'apo.12': '4jhw_tail_f', 'apo.13': '5udc_monomer_f', 
    'apo.14': '5udc_tail_monomer_f', 'apo.15': '5udc_f', 'apo.16': '5udc_tail_f'}

infile_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_output/"
outfile_prefix = "/data/chodera/zhangi/vir_collaboration/data/em_figures/"

for k, v in sorted(d.items()):
    print("Trying: ", v)
    state_data = infile_prefix + k + '.50ns.csv'
    ref_pdb = infile_prefix + k + '.50ns.minimized.pdb'
    traj = infile_prefix + k + '.50ns.solute.dcd'
    if 'holo' in k: # For restarted sims only, remove this in future analysis
        state_data = infile_prefix + k + '.50ns.combined.csv'
        ref_pdb = infile_prefix + k + '.50ns.minimized.pdb'
        traj = infile_prefix + k + '.50ns.combined.solute.dcd'
    plot_all(state_data, v)
    plot_rmsd(traj, ref_pdb, v)
    print('done')