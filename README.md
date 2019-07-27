# antibody-mutations
Assessing the impact of mutations on antibody-antigen binding affinity

## Clean PDBs
### Scripts
clean_pdb.ipynb, edit_seqres.ipynb
### Description of clean_pdb.ipynb
Use this jupyter notebook to clean up PDBs.
1. Preprocess PDBs -- e.g. If there is a loop you wish to remove, manually remove the loops from ATOM lines and MISSING RESIDUES. Use edit_seqres.ipynb to programttically remove the loop from the SEQRES.
* Description of edit_seqres.ipynb: This script calls 3 functions to manipulate a SEQRES: 1) delete lines 2) fill SEQRES lines that are not full (13 residues)3) reindex SEQRES lines so that indexes are consecutive. Output: PDB file with desired SEQRES sequence removed
2. Use MDTraj to renumber residues in each chain such that there are no residue numbers with letters (due to insertion codes) and there are gaps in residues that corresponding to the missing residues we want added back in.
3. Use PDBFixer to add missing residues and remove heterogens. Note: If the missing residues were part of a long terminal fragment (i.e. > 10 residues), then they will not be added them back in.
### Input files
data/4jhw/4jhw.pdb and data/5udc/5udc.pdb
### Output files
* data/5udc/5udc_noloop_noseqres.pdb — manually removed loop from ATOM lines and MISSING RESIDUES and programmatically removed loop from SEQRES.
* data/5udc/5udc_mdtraj_noloop_noseqres_nogap.pdb — mdtraj intermediate (with residue numbers renumbered, gaps edited, and SEQRES copied from 5udc_noloop_noseqres.pdb)
* data/5udc/5udc_clean_nolongterms_noseqres_nogap.pdb — no long (> 10 residues) terminal fragments added back in. Ran PDBFixer on 5udc_mdtraj_noloop_noseqres_nogap.pdb
* data/4jhw/4jhw.pdb — raw PDB
* data/4jhw/4jhw_mdtraj.pdb — mdtraj intermediate (with residue numbers renumbered, gaps edited, and SEQRES copied from 4jhw.pdb)
* data/4jhw/4jhw_clean.pdb — no long (> 10 residues) terminal fragments added back in
