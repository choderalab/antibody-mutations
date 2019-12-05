
## PREP STRUCTURE/TOPOLOGY
# Create topology
gmx pdb2gmx -f 5udc_final_v2_refmac1_clean_tail_variable_monomer.pdb -o conf.pdb -ff amber99sb-ildn -water tip3p

## SET UP SIMULATION
# Define simulation box
gmx editconf -f conf.pdb -o box.pdb -bt dodecahedron -d 0.9

# Solvate simulation box
gmx solvate -cp box.pdb -cs spc216 -o water.pdb -p topol.top

# Assemble atomic level description of system (.tpr file) 
gmx grompp -f 01_ions.mdp -c water.pdb -p topol.top -o ions.tpr

# Add chloride atoms
gmx genion -s ions.tpr -conc 0.150 -neutral -o ions.pdb -nname CL -pname NA -p topol.top

# Assemble .tpr file again (this time for input to mdrun)
gmx grompp -f 02_em.mdp -c ions.pdb -p topol.top

## RUN SIMULATION
# Run energy minimization
gmx mdrun -v -c em.pdb

# Equilibrate the system 
echo 0 | gmx trjconv -s topol.tpr -f em.pdb -o em.pdb -ur compact -pbc mol
gmx grompp -f 03_eq.mdp -c em.pdb -p topol.top
gmx mdrun -v -ntomp 2 -nice 0 -c eq.gro

# Run MD
gmx grompp -f 04_md.mdp -c eq.gro -p topol.top
gmx mdrun -v -ntomp 2 -nice 0 -c md.gro
