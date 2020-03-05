# Nanodisc-Builder
Protocol for constructing nanodiscs in silico

Usage:

> ./Build_MSP_monomer.py -f base.fasta -o base

> ./Assemble_two_monomers.py -m base.pdb -f base.fasta -w 

This construct an atomistic model of two MSPs as a dimer, where lipids and/or a membrane protein can be inserted.

The COMMANDS_default and COMMANDS_Uni are both 'wrapper' bash scripts, which construct the nanodisc dimer using the python scripts and add lipids and solvate using insane. 

The COMMANDS_default uses the deafult settings:
Interface as Left/Left 5/5, &
The distance between monomers as 12 Å

Within the two 'wrapper' scripts settings are set in top. 
Both scripts are for using the Martini and Charmm36 FF with the GROMACS simulation software. 

For help or documentation of scripts run:

> ./Build_MSP_monomer.py -h
> ./Assemble_two_monomers.py -h

> ./COMMANDS_Uni --help
> ./COMMANDS_default --help
