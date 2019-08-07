# antibody-mutations
Assess the impact of mutations on antibody-antigen binding affinity

## Clean PDBs
### Scripts
clean_pdb.ipynb, edit_seqres.ipynb, split_chain_and_cap.ipynb
### Description of protocol
1. Preprocess PDBs -- Use edit_seqres.ipynb to programatically split delete lines and/or split chains in the SEQRES.
    * Description of edit_seqres.ipynb: This script calls 3 functions to manipulate a SEQRES: 1) delete lines (or split_chains) 2) fill SEQRES lines that are not full (13 residues)3) reindex SEQRES lines so that indexes are consecutive. Output: PDB file with desired SEQRES sequence removed
2. Split and cap F protein chain -- Use split_chain_and_cap.ipynb
    * Description of split_chain_and_cap.ipynb: This script uses the OpenMM Topology object to create a new topology with the F protein chain(s) split into two (because it contains a long missing loop that was not added back in) and caps them using PDBFixer.
3. Manually copy the edited SEQRES to the PDB with the F protein split and capped.
4. Clean PDB -- Use clean_pdb.ipynb
    * Description of clean_pdb.ipynb: This script uses PDBFixer to add missing residues and remove heterogens. Note: If the missing residues were part of a long terminal fragment (i.e. > 10 residues), then they will not be added them back in.

### Input files
data/4jhw/4jhw.pdb and data/5udc/5udc.pdb
### Output files
Note: the files for each structure are in order from input to output and accumulate changes, e.g. 4jhw_noloop_noseqres.pdb contains changes from 4jhw_noloop.pdb)

5udc
* data/5udc/5udc.pdb - raw PDB
* data/5udc/5udc_noloop.pdb - deleted residues that are not present in 4jhw, 4jha, 5k6f
* data/5udc/5udc_noloop_seqressplit.pdb - split chain in SEQRES near missing loop that will not be added back in
* data/5udc/5udc_splitchain.pdb - F protein chain split into two chains
* data/5udc/5udc_splitchain_capped.pdb - chains capped using PDBFixer
* data/5udc/5udc_clean.pdb — missing residues not in long (> 10 residues) terminal fragments added back in

4jhw
* data/4jhw/4jhw.pdb — raw PDB
* data/4jhw/4jhw_noloop.pdb - deleted residues that are not present in 5udc, 4jha, 5k6f
* data/4jhw/4jhw_noloop_seqressplit.pdb - split chain in SEQRES near missing loop that will not be added back in
* data/4jhw/4jhw_splitchain.pdb - F protein chain split into two chains
* data/4jhw/4jhw_splitchain_capped.pdb - chains capped using PDBFixer
* data/4jhw/4jhw_clean.pdb — missing residues not in long (> 10 residues) terminal fragments added back in
### Example protocol
4jhw
* Load 4jhw.pdb into edit_seqres.ipynb to delete residues and split chain in the SEQRES near missing loop that will not be added back in. Write this to 4jhw_noloop.pdb (after deleting residues), then 4jhw_noloop_seqressplit.pdb.
* Load 4jhw.pdb into split_chain_and_cap.ipynb. Save topology and positions with chain split to 4jhw_splitchain.pdb.
* Load 4jhw_splitchain.pdb into split_chain_and_cap.ipynb. Follow PDBFixer protocol (findMissingResidues(), findMissingAtoms(), addMissingAtoms()) to add missing terminal atoms. Save topology and positions with chains capped to 4jhw_splitchain_capped.pdb.
* Manually copy the SEQRES from 4jhw_seqressplit.pdb to the top of 4jhw_splitchain_capped.pdb.
* Load 4jhw_splitchain_capped.pdb into clean_pdb.ipynb. Follow PDBFixer protocol (findMissingResidues(), remove terminal missing fragments from fixer.missingResidues, findNonStandardResidues(), removeHeterogens, findMissingAtoms(), addMissingAtoms()). Write to 4jhw_clean.pdb.
