# antibody-mutations
Assessing the impact of mutations on antibody-antigen binding affinity

## Clean PDBs
### Description
Use this jupyter notebook to clean up PDBs.
1. Use MDTraj to renumber residues in each chain such that there are no residue numbers with letters (due to insertion codes).
2. Use PDBFixer to add missing residues and remove heterogens. Note: If the missing residues were part of a long terminal fragment (i.e. > 10 residues), then they will not be added them back in.
Note: For 4jhw, the F protein SEQRES was missing a segment of residues that was present in 5udc, so I manually edited the SEQRES and MISSING RESIDUES section to include this segment.
### Script
clean_pdb.ipynb
### Input files
data/4jhw/4jhw.pdb and data/5udc/5udc.pdb
### Output files
* data/5udc/5udc_mdtraj.pdb — mdtraj intermediate (with residue numbers renumbered and SEQRES copied from 5udc.pdb)
* data/5udc/5udc_clean_noterm.pdb — has no terminal fragments added in at all, not even when its just a singular missing residue… is this what we want?
* data/5udc/5udc_clean_nolongterms.pdb — only short terminal fragments added back in
* data/4jhw/4jhw_add_missing.pdb — SEQRES and missing residues sections altered to include middle missing fragment residues that are present in 5udc but not 4jhw
* data/4jhw/4jhw_clean.pdb — no long terminal fragments added back in
* data/4jhw/4jhw_mdtraj_add_missing.pdb — mdtraj intermediate (with residue numbers renumbered and SEQRES copied from 4jhw.pdb), with missing residues (from middle missing fragment) present in 5udc but not 4jhw added in 
