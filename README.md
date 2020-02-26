# antibody-mutations
Assess the impact of mutations on antibody-antigen binding affinity

### Software package versions
OpenMM: 7.3.1

OpenMMTools: 0.16.1

CUDA: 8.0

## Clean PDBs
### Scripts
clean_pdb_rerefined.ipynb, clean_pdb_4jha.ipynb, edit_seqres.ipynb, split_chain_and_cap.ipynb
### Instructions (for 5udc and 4jhw)
1. Preprocess SEQRES -- Use edit_seqres.ipynb to programatically split delete lines and/or split chains in the SEQRES.
    * Description of edit_seqres.ipynb: This script calls 3 functions to manipulate a SEQRES: 1) delete lines (or split_chains) 2) fill SEQRES lines that are not full (13 residues) 3) reindex SEQRES lines so that indexes are consecutive. Output: PDB file with desired SEQRES sequence removed
2. Split F protein chain -- Use split_chain.ipynb
    * Description of split_chain.ipynb: This script uses the OpenMM Topology object to create a new topology with the F protein chain(s) split into two (because it contains a long missing loop that was not added back in). Terminal atoms added to ends of chains using PDBFixer.
3. Manually copy CRYST1 line and the edited SEQRES to the PDB with the F protein split.
4. Clean PDB -- Use clean_pdb.ipynb
    * Description of clean_pdb.ipynb: This script uses PDBFixer to add missing residues and remove heterogens. Note: If the missing residues were part of a long terminal fragment (i.e. > 10 residues), then they will not be added them back in.
5. Add CRYST1 line from the original PDB if the final PDB doesn't contain it.

### Input files
data/4jhw/4jhw.pdb and data/5udc/5udc.pdb
### Output files
Note: the files for each structure are in order from input to output and accumulate changes, e.g. 4jhw_noloop_noseqres.pdb contains changes from 4jhw_noloop.pdb)

5udc
* data/5udc/5udc.pdb - raw PDB
* data/5udc/5udc_noloop.pdb - deleted residues that are not present in 4jhw, 4jha, 5k6f
* data/5udc/5udc_noloop_seqressplit.pdb - split chain in SEQRES near missing loop that will not be added back in
* data/5udc/5udc_splitchain.pdb - F protein chain split into two chains
* data/5udc/5udc_splitchain_capped.pdb - terminal atoms added to ends of chains using PDBFixer
* data/5udc/5udc_clean.pdb — missing residues not in long (> 10 residues) terminal fragments added back in

4jhw
* data/4jhw/4jhw.pdb — raw PDB
* data/4jhw/4jhw_noloop.pdb - deleted residues that are not present in 5udc, 4jha, 5k6f
* data/4jhw/4jhw_noloop_seqressplit.pdb - split chain in SEQRES near missing loop that will not be added back in
* data/4jhw/4jhw_splitchain.pdb - F protein chain split into two chains
* data/4jhw/4jhw_splitchain_capped.pdb - terminal atoms added to ends of chains using PDBFixer
* data/4jhw/4jhw_clean.pdb — missing residues not in long (> 10 residues) terminal fragments added back in
### Example
4jhw
* Load 4jhw.pdb into edit_seqres.ipynb to delete residues and split chain in the SEQRES near missing loop that will not be added back in. Write this to 4jhw_noloop.pdb (after deleting residues), then 4jhw_noloop_seqressplit.pdb.
* Load 4jhw.pdb into split_chain.ipynb. Use the split_chain function to save topology and positions with chain split to 4jhw_splitchain.pdb. Then use the cap_chain function to save topology and positions with missing terminal atoms added to 4jhw_splitchain_capped.pdb.
* Manually copy the SEQRES from 4jhw_noloop_seqressplit.pdb to the top of 4jhw_splitchain_capped.pdb.
* Load 4jhw_splitchain_capped.pdb into clean_pdb.ipynb. Follow PDBFixer protocol (findMissingResidues(), remove terminal missing fragments from fixer.missingResidues, findNonStandardResidues(), removeHeterogens, findMissingAtoms(), addMissingAtoms()). Write to 4jhw_clean.pdb.
### Instructions for 5k6f
For this structure, remove the 99-109 region (and the linker) and then split the chain at this region, using split_chain_and_cap.ipynb. There are no missing residues so no need to clean with PDBFixer.
### Instructions 4jha
For this structure, add missing residues based on 4jhw antibody using clean_4jha.ipynb. There is no need to split chains here.
### Instructions for 5udc trimer F protein + single antibody
1. Determine antibody with best electron density and use UCSF Chimera to remove chains from the unwanted antibodies.
2. Add CRYST1 line from the original PDB if the final PDB doesn't contain it
### Instructions for 4jhw trimer F protein + single antibody
1. Copy 4jhw_final_v2_refmac1_clean.pdb to new file 4jhw_final_v2_refmac1_clean_trimer.pdb.
2. Load 4jhw_final_v2_refmac1_clean_trimer.pdb as "Biological Assembly" in UCSF Chimera
* Render biological assembly: Model panel > biological unit
* Combine the assembly into one model: Copy/combine
* Save the PDB: File > Save PDB 
3. Remove chains of unwanted antibodies in Chimera
* Change chain ids to match 5udc trimer F protein chains: Tools > Structure editing > Change chain IDs
* Save the PDB: File > Save PDB (renumbered/holo/nonoverlay/4jhw_final_v2_refmac1_clean_trimer_single_ab.pdb)
4. Add CRYST1 line from the original PDB if the final PDB doesn't contain it

## Remove F protein head or Antibody constant
### Scripts
delete_residues.ipynb
### Description
This script deletes residues at the n-terminus or c-terminus as specified by the user. It caps the chain at the terminus with deleted residues.

## Overlay high-resolution structures onto low-resolution structures
0. Load PDBs into UCSF Chimera
   * File > Open > Select PDB file
1. Overlay
   * Tools > Structure comparison > Matchmaker (Use default settings)
2. Delete low res reference chains
   * Select > Chain 
   * Actions > Atoms/Bonds > Delete 
3. Merge models
   * Tools > General controls > Model panel > Select the two models, choose “copy/combine"
4. Change Chain IDs so that they correspond to chain IDs in the low res holo structure
   * Tools > Structure editing > Change chain IDs
5. File > Save PDB
   * Save with respect to low res structure

## Extract monomer from 5udc
0. Load PDB into UCSF Chimera
   * File > Open > Select PDB file
1. Delete all chains for two of the monomers
   * Select > Chain 
   * Actions > Atoms/Bonds > Delete 
   * Note: To select multiple chains at a time (e.g. chains F through L), enter into Chimera's command line: `select  :.f-l`
   * Note: I kept H, L, F chains in the original PDB (aka chains A-D in the cleaned PDB) because these required adding the least number of missing loop residues
2. File > Save PDB

Note: To extract the F protein or antibody, follow this same procedure, but select the relevant chains.

## Run energy minimization simulations
### Scripts
run_em.py, run_em.sh, remove_solvent_from_traj.py, join_restarted.py
### Instructions
(to be run on lilac)
1. Run `module load cuda/9.2`
2. Activate conda environment
3. Navigate to the directory containing run_em.py and run_em.sh and run `bsub < run_em.sh`
4. After jobs have completed, remove solvent from trajectories using remove_solvent_from_traj.py
5. If there were jobs that were restarted, join the original and restarted trajectories and CSVs using join_restarted.py

## Restart a simulation
### Scripts
restart_em.py, restart_em.sh
### Instructions
1. Run `module load cuda/9.2`
2. Activate conda environment
3. Navigate to the directory containing restart_em.py and restart_em.sh and run `bsub < restart_em.sh`
4. After jobs have completed, remove solvent from trajectories using remove_solvent_from_traj.py

## Analyze simulations
### Scripts
analyze_em.py, analyze_em_interface_residues.py, analyze_em_interface_rmsd.py
### Instructions
1. Run analyze_em.py to create plots of state data (PE, KE, T, V) and RMSD (after aligning on whole molecule).
2. Run analyze_em_interface_residues.py to create dataframes of interface residues in each chain.
3. Run analyze_em_interface_rmsd.py to plot RMSDs (interface, F protein, antibody) after aligning different components (interface, F protein, antibody).

## Saving simulations as a movie:
(In PyMOL)

Rafal’s tips:
1. Display > Quality > Maximum performance
2. Display > Color Space > CMYK
3. set ray_trace_mode, 3
4. show spheres
5. set sphere_transparency, 0.85
6. set movie_fps, 60
7. smooth
8. zoom
9. Center it on the right view of the interface
10. Smooth (can do this multiple times, will remove the faster movements more and more)
11. Save as pymol sessions (if need to align complexes in different movies, can just load other complex into this session and then superpose/align the other complex then save as movie)
12. File > Export as movie.. MPEG (make the window as big of a square as possible —> 1600x1600 for high quality)
Save every 0.1-0.5 ns
