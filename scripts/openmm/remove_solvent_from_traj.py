import mdtraj as md
import os

directory = "/data/chodera/zhangi/vir_collaboration/data/em_output/"
for filename in os.listdir(directory):
	if ".50ns.dcd" in filename:
		topology_file = directory + filename[:-3] + "minimized.pdb"
		output_file = directory + filename[:-3] + "solute.dcd"
		traj = md.load(directory + filename, top= topology_file)
		traj_r = traj.remove_solvent()
		traj_r.save_dcd(output_file)