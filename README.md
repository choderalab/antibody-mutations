# antibody-mutations
Assess the impact of mutations on antibody-antigen binding affinity

## Clean PDBs
### Scripts
clean_pdb.ipynb, edit_seqres.ipynb, split_chain_and_cap.ipynb
### Description of protocol
1. Preprocess PDBs -- e.g. If there is a loop you wish to remove, manually remove the loops from ATOM lines and MISSING RESIDUES. Use edit_seqres.ipynb to programatically remove the loop from the SEQRES.
    * Description of edit_seqres.ipynb: This script calls 3 functions to manipulate a SEQRES: 1) delete lines 2) fill SEQRES lines that are not full (13 residues)3) reindex SEQRES lines so that indexes are consecutive. Output: PDB file with desired SEQRES sequence removed
2. Clean PDB -- Use clean_pdb.ipynb
    * Description of clean_pdb.ipynb: This script uses MDTraj to renumber residues in each chain such that there are no residue numbers with letters (due to insertion codes) and there are gaps in residues that corresponding to the missing residues we want added back in. Then, it uses PDBFixer to add missing residues and remove heterogens. Note: If the missing residues were part of a long terminal fragment (i.e. > 10 residues), then they will not be added them back in.
3. Split and cap F protein chain -- Use split_chain_and_cap.ipynb
    * Description of split_chain_and_cap.ipynb: This script uses the OpenMM Topology object to create a new topology with the F protein chain(s) split into two (because it contains a long missing loop that was not added back in) and caps them using PDBFixer.
### Input files
data/4jhw/4jhw.pdb and data/5udc/5udc.pdb
### Output files
Note: the files for each structure are in order from input to output and accumulate changes, e.g. 4jhw_noloop_noseqres.pdb contains changes from 4jhw_noloop.pdb)

5udc

* data/5udc/5udc_noloop_noseqres.pdb — manually removed loop from ATOM lines and MISSING RESIDUES and programmatically removed loop from SEQRES.
* data/5udc/5udc_mdtraj_noloop_noseqres_nogap.pdb — mdtraj intermediate (with residue numbers renumbered, gaps edited, and SEQRES copied from 5udc_noloop_noseqres.pdb)
* data/5udc/5udc_clean_nolongterms_noseqres_nogap.pdb — no long (> 10 residues) terminal fragments added back in. Ran PDBFixer on 5udc_mdtraj_noloop_noseqres_nogap.pdb
* data/5udc/5udc_clean_nolongterms_noseqres_nogap.pdb — no long (> 10 residues) terminal fragments added back in. Ran PDBFixer on 5udc_mdtraj_noloop_noseqres_nogap.pdb


4jhw
* data/4jhw/4jhw.pdb — raw PDB
* data/4jhw/4jhw_noloop.pdb - missing loops that should not be added back in manually removed from MISSING RESIDUES
* data/4jhw/4jhw_noloop_noseqres.pdb - missing loops that should not be added back in removed from SEQRES
* data/4jhw/4jhw_mdtraj.pdb — mdtraj intermediate (with residue numbers renumbered, gaps edited, and SEQRES copied from 4jhw.pdb)
* data/4jhw/4jhw_clean.pdb — missing residues not in long (> 10 residues) terminal fragments added back in
* data/4jhw/4jhw_clean_splitchain.pdb - F protein chain split into two chains
* data/4jhw/4jhw_clean_splitchain_capped.pdb - chains capped using PDBFixer
### Example protocol
4jhw
* Make a duplicate copy of the raw pdb (4jhw_noloop.pdb) and manually remove the loop(s) that will not be added back in from the MISSING RESIDUES section
* Load 4jhw_noloop.pdb into edit_seqres.ipynb to remove the loop(s) that will not be added back in from the SEQRES. Write this to 4jhw_noloop_noseqres.pdb.
* Load raw pdb into clean_pdb.ipynb. Use mdtraj to renumber residues, maintaining gaps for missing residues. Write topology and trajectory to 4jhw_mdtraj.pdb.
* Manually copy the SEQRES from 4jhw_noloop_noseqres.pdb to the top of 4jhw_mdtraj.pdb.
* Load 4jhw_mdtraj.pdb into clean_pdb.ipynb. Follow PDBFixer protocol (findMissingResidues(), remove terminal missing fragments from fixer.missingResidues, findNonStandardResidues(), removeHeterogens, findMissingAtoms(), addMissingAtoms()). Write to 4jhw_clean.pdb.
* Load 4jhw_clean.pdb into split_chain_and_cap.ipynb. Save topology and positions with chain split to 4jhw_clean_splitchain.pdb.
* Load 4jhw_clean_splitchain.pdb into split_chain_and_cap.ipynb. Follow PDBFixer protocol (findMissingResidues(), findMissingAtoms(), addMissingAtoms()) to add missing terminal atoms. Save topology and positions with chains capped to 4jhw_clean_splitchain_capped.pdb.
