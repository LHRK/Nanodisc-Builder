#!/bin/bash

### Bash script for MIN + EG of Coarse Grained ND systems ###


gmx grompp -f min.mdp -c ${base}_cg-solvent.gro  -p ${base}_cg-solvent.top -n ${base}_cg-solvent.ndx -maxwarn 1 -o ${base}_cg-solvent-min.tpr
gmx mdrun -v -deffnm ${base}_cg-solvent-min

gmx grompp -f eq-NVT.mdp -c ${base}_cg-solvent-min.gro -p ${base}_cg-solvent.top -n ${base}_cg-solvent.ndx -maxwarn 1 -o ${base}_cg-solvent-eq-NVT.tpr
gmx mdrun -v -deffnm ${base}_cg-solvent-eq-NVT

gmx grompp -f eq-NPT-0.005.mdp -c ${base}_cg-solvent-eq-NVT.gro -p ${base}_cg-solvent.top -n ${base}_cg-solvent.ndx -maxwarn 1 -o ${base}_cg-solvent-eq-NPT-0.005.tpr
gmx mdrun -v -deffnm ${base}_cg-solvent-eq-NPT-0.005

gmx grompp -f eq-NPT-0.01.mdp -c ${base}_cg-solvent-eq-NPT-0.005.gro -p ${base}_cg-solvent.top -n ${base}_cg-solvent.ndx -maxwarn 1 -o ${base}_cg-solvent-eq-NPT-0.01.tpr
gmx mdrun -v -deffnm ${base}_cg-solvent-eq-NPT-0.01 

gmx grompp -f eq-NPT-0.020.mdp -c ${base}_cg-solvent-eq-NPT-0.01.gro -p ${base}_cg-solvent.top -n ${base}_cg-solvent.ndx -maxwarn 1 -o ${base}_cg-solvent-eq-NPT-0.020.tpr
gmx mdrun -v -deffnm ${base}_cg-solvent-eq-NPT-0.020 

gmx grompp -f md.mdp -c  ${base}_cg-solvent-eq-NPT-0.020.gro -p ${base}_cg-solvent.top -n ${base}_cg-solvent.ndx -maxwarn 1 -o ${base}_cg-solvent-md.tpr

