## PREP STRUCTURE/TOPOLOGY
# Create morph/hybrid structure
# python ../../pmx/scripts/mutate.py -f 5udc_final_v2_refmac1_clean_tail_variable_monomer_CH3.pdb -o mut.pdb -ff amber99sbmut -script test_mutation.txt

# Create topology
gmx pdb2gmx -f 5udc_final_v2_refmac1_clean_tail_variable_monomer_CH3.pdb -ff amber99sb-star-ildn-mut -water tip3p -heavyh -norenum -o conf.pdb 
# gmx pdb2gmx -f mut.pdb -ff amber99sb-star-ildn-mut -water tip3p -heavyh -norenum -o mut_pdb2gmx.pdb 

# Add morphing paramters to topology
# python ../../pmx/scripts/generate_hybrid_topology.py -p topol.top -o hybrid.top -ff amber99sbmut

## SET UP SIMULATION
# Define simulation box
gmx editconf -f conf.pdb -o box.pdb -bt cubic -d 1
# gmx editconf -f mut_pdb2gmx.pdb -o box.pdb -bt cubic -d 1

# Solvate simulation box
gmx solvate -cp box.pdb -cs spc216 -o water.pdb -p topol.top
# gmx solvate -cp box.pdb -cs spc216 -o water.pdb -p hybrid.top

# Assemble atomic level description of system (.tpr file) 
gmx grompp -f 00_ions.mdp -c water.pdb -p topol.top -norenum -o ions.tpr 
# gmx grompp -f 00_ions.mdp -c water.pdb -p hybrid.top -norenum -o ions.tpr 

# Add chloride atoms
gmx genion -s ions.tpr -conc 0.150 -neutral -o ions.pdb -nname CL -pname NA -p topol.top
# gmx genion -s ions.tpr -conc 0.150 -neutral -o ions.pdb -nname CL -pname NA -p hybrid.top

# Assemble .tpr file again (this time for input to mdrun)
gmx grompp -f 01_eminim.mdp -c ions.pdb -p topol.top -norenum
# gmx grompp -f 01_eminim.mdp -c ions.pdb -p hybrid.top -norenum

## RUN SIMULATION
# Run energy minimization
# gmx mdrun -v -c em.pdb -x em.xtc -nt 4
gmx mdrun -v -c em.pdb -x em.xtc -cpo -nt 4 
# gmx mdrun -v -c em.pdb -x em.xtc -cpi state.cpt -append -nt 4 # Restart

# Equilibrate the system 
echo 0 | gmx trjconv -s topol.tpr -f em.pdb -o em.pdb -ur compact -pbc mol
gmx grompp -f 02_equilibration.mdp -c em.pdb -r em.pdb -p topol.top -norenum
# gmx grompp -f 02_equilibration.mdp -c em.pdb -r em.pdb -p hybrid.top -norenum
gmx mdrun -v -c eq.pdb -o eq.trr -x eq.xtc -cpo eq.cpt -e eq.edr -g eq.log -nt 4

# Run MD at equilibrium
gmx grompp -f 03_equilibrium_prod.mdp -c eq.pdb -r eq.pdb -p topol.top -norenum
# gmx grompp -f 03_equilibrium_prod.mdp -c eq.pdb -p hybrid.top -norenum
gmx mdrun -v -c md_e.pdb -o md_e.trr -x md_e.xtc -cpo md_e.cpt -e md_e.edr -g md_e.log -nt 4 # do i need to specify xvg files?
