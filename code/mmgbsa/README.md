# antibody-mutations
Assess the impact of mutations on antibody-antigen binding affinity

### File descriptions
* mmgbsa.py
    * Script for running MM/GBSA. (See instructions.txt)
* input.txt
    * input_file_PL_WT
        * Path to input file containing WT protein-ligand complex on which to compute MM/GBSA
    * input_file_PL_mutant
        * Path to input file containing mutant protein-ligand complex on which to compute MM/GBSA
    * add_hydrogens
        * Specify True if you want to add hydrogens to each structure.
        * Specify False if the structure already contains hydrogens.
    * output_directory
        * Path to output directory
    * name
        * Name of complex (this will be used in output file names)
    * protein_chains
        * Comma-separated chain ids of the protein chains in the complex (as named in input_file)
    * ligand_chains
        * Comma-separated chain ids of the ligand chains in the complex (as named in input_file)
    * iterations
        * Number of MM/GBSA estimates per complex.
* WT.pdb
    * 6nb8_2ghv_sars structure (un-mutated, capped)


### Use
1. Modify parameters in `input.txt` or create new .txt file with desired parameters.
2. Run the following in command line: `python mmgbsa.py input.txt`