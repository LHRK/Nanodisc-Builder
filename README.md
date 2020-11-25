# Nanodisc-Builder
Protocol for constructing nanodiscs in silico

## Installation 
It is recommended to make a virtual environment with eg Anaconda like so:  

> conda create -n ND_Builder python=2.7\
> conda activate ND_Builder\
> conda install numpy MDAnalysis rmsd Biopyhton argparse re os\
> pip install PeptideBuilder\
> pip install insane\
> conda deactivate

### Requirements
Oterwise simply make sure to have installed the following modules:\
- numpy
- argsparse
- os
- re
- MDAnalysis
- PeptideBuilder
- Biopython
- rmsd 


and possibly insane as well:
> pip install insane\
Insane is however, only compatible with python 2.7, as is the martinize script as well. 


# Python version
The python script:
- Build_MSP_monomer.py
- Assemble_two_monomers.py
- Calc_number_lipids.py
are all compatible with python2.7 and python3

The Insane and Martinize scripts are however only compatible in python2.7. 
The Martinize script can also be found at the github page: https://github.com/cgmartini/martinize.py, where a python 3 version is avaliable upon request. 

For all the python scripts, additional flags and documentation can be viewed upon running the script with the flag -h. 

# Usage

## Constructing the atomistic MSP dimer 
Frist go into the generated virtual eenvironment.
> conda activate ND_Builder

Use the Build_MSP_monomer.py to construc the MSP monomer given the fasta sequence.\
The output structure may appear to have clashed at the kinks between the tandem repeats. This is however solved upon a minimization.

> ./Build_MSP_monomer.py -f base.fasta -o base

For assembling the MSP dimer, use the Assemble_two_monomers.py script

> ./Assemble_two_monomers.py -m base.pdb -f base.fasta -w

For using default setting use the -w flag. This includes assembling the nanodisc with the LL5/5 registry, given the helix5 sequence is found in the input MSP sequence.\
A distance of 12 Å between the MSP monomers is the default setting.

For changing the registry the following flag can be used:\
--if : for selecting either LL (left-left) or RR (right-right) registry\
-t   : Changes the distance between the two MSP monomers i Angstrom. Default is 12 Å.

For changing the tandem or sequence to superimpose on the following flags can be used:\
--r1\
--r2\
-H

The -H flag can be used to select a tandem repeat for assembling the MSP dimer.\
Eg. -H 'H1'\
The avaliabel tandem repeats are listed below.\
If the tandem repeat is not in the list, the flags --r1 and --r2 can be used to state the beginning and end resid of the desired tandem repeat.\
Eg. --r1 105 --r2 115

#### Tandem repeats
H1 : LKLLDNWDSVTSTFSKLREQLG\
H2 : PVTQEFWDNLEKETEGLRQEMS\
H3 : KDLEEVKAKVQ\
H4 : PYLDDFQKKWQEEMELYRQKVE\
H5 : PLRAELQEGARQKLHELQEKLS\
H6 : PLGEEMRDRARAHVDALRTHLA\
H7 : PYSDELRQRLAARLEALKENGG\
H8 : ARLAEYHAKATEHLSTLSEKAK\
H9 : PALEDLRQGLL\
H10 : PVLESFKVSFLSALEEYTKKLNTQ


## Coarse grain and insert lipid and possible a membrane protein
The script martinize and insane can then following be used for coarse graining the MSP dimer and inserting lipids and possible a membrane protein.

For insane an estimated radius of the ND is needed.\
This can be calculated with the python script
> ./Calc_number_lipids.py -h\
Eg\
> ./Calc_number_lipids.py -f base.fasta -a 70 -p 0  


The -a flag is the APL estimate for the lipid type in Angstrom. This is only for single lipid mixtures.\
The -p flag is the approximate transmembrane area, which will be occupied by a possible membrane protein in the ND.  

The script output a Suggestions.txt file with an estimated radius of the disc, along with estimated number of lipids er leaflet, based on the given APL.

### Mixed lipid compositions and embedding of membrane protein 
Insane can be used for mixed lipid compositions in the ND, but with and without a possible embedded membrane protein. 

### Circularized NDs
The script Corr_itp_circularized.sh is for correcting the itp file for the coarse grained MSP dimer, making sure the terminals are covalently linked together.\
The format is for the GROMACS software. 

## Minimization and equlibration
The Run.sh script can be used for suggested steps of minimization and equlibration.\
For GROMACS software version 2018.2 and up the Run_ver18.sh script is valid.  

## Backmapping to atomistic
For converting the coarse grained ND to atomistic scale the backward.py script can be used.\
It is however, recommended not to backmap the waters and ions, but only the ND complex and then re-solvate and neutralize.\
For the GROMACS software the Backmap_re-solvate.sh script can be used.\
Edit in the framed box within the script.

# Wrapper scripts for automazation 
The COMMANDS_default and COMMANDS_Uni are both 'wrapper' bash scripts, which construct the nanodisc dimer using the python scripts and add lipids and solvate using insane.

The COMMANDS_default uses the deafult settings:\
Interface as Left/Left 5/5, &\
The distance between monomers as 12 Å

Within the two 'wrapper' scripts settings are set in top.\
Both scripts are for using the Martini and Charmm36 FF with the GROMACS simulation software.\

> ./COMMANDS_Uni --help\
> ./COMMANDS_default --help

For using the wrapper scripts with more complex lipid compositions used the COMMANDS_Uni script and edit the option:\
lipid_type='DLPC' in the script to eg. lipid_type='-l POPC:90 -l POPG:10 -u POPC:90 -l POPG:10'   
This line will be given directly to insane, which will then insert 90% POPC and 10% POPG in both the lower (-l) and upper (-u) leaflet.\
Hence the two leaflets can be constructed differently as well, by using the -l and -u flags provided to the insane script.

# Examples

## A simple empty ND

## An empty ND with complex lipid composition

## A ND with a GPCR embedded - simple lipid composition

## A ND with a GPCR embedded - complex lipid composition 
