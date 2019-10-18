import mdtraj as md
import pandas as pd

# Specify file names that need to be joined
file_prefix = '/data/chodera/zhangi/vir_collaboration/data/em_output/'
files = ['holo.nonoverlay.5', 'holo.nonoverlay.6', 'holo.overlay.4jha.1', 'holo.overlay.4jha.5', 'holo.overlay.5k6f.1',
'holo.overlay.5k6f.5', 'holo.overlay.5k6f.9', 'holo.overlay.5k6f.10']

# Choose file paths for original and restarted trajectories and csvs
for file in files:
	print("Trying ", file)
	original_traj_file = file_prefix + file + '.50ns.solute.dcd'
	original_ref_file = file_prefix + file + '.50ns.minimized.pdb'
	restarted_traj_file = file_prefix + file + '.50ns.restarted.solute.dcd'
	original_data_file = file_prefix + file + '.50ns.csv'
	restarted_data_file = file_prefix + file + '.50ns.restarted.csv'

	## JOIN CSVs
	print("\tJoining CSVs")
	df_orig = pd.read_csv(original_data_file)
	df_re = pd.read_csv(restarted_data_file)
	n_orig_rows = df_orig.shape[0]
	df_orig = df_orig.append(df_re.iloc[:1000-n_orig_rows])
	df_orig.to_csv(file_prefix + file + '.50ns.combined.csv', index=False)

	### JOIN TRAJECTORIES
	print("\tLoading topology PDB")
	# Load topology PDB
	original_ref = md.load_pdb(original_ref_file)
	print("\tRemoving solvent from topology PDB")
	original_ref_solute = original_ref.remove_solvent()

	# Load original trajectory
	print("\tLoading original trajectory")
	original_traj = md.load_dcd(original_traj_file, top=original_ref_solute)

	# Load restarted trajectory, slicing to keep only the number of additional frames needed
	print("\tLoading restarted trajectory")
	restarted_traj = md.load(restarted_traj_file, top=original_ref_solute)
	restarted_traj = restarted_traj[:1000-original_traj.n_frames]

	# Join the trajectories
	print("\tJoining trajectories and saving them")
	joined_traj = md.join([original_traj, restarted_traj])
	joined_traj.save(file_prefix + file + ".50ns.combined.solute.dcd")