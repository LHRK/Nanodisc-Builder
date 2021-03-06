# Generic protocol for constructing Nanodiscs in silico
The protocol introduces a novel approach for constructing MSPs in silico.
Integrated with Insane and Martinize scripts, the protocol can be used to easily construct nanodisc systems.

## Installation 
It is recommended to make a virtual environment with eg Anaconda like so:  
(Supportive python verions: 3)
> conda create -n ND_Builder\
> conda activate ND_Builder\
> conda install numpy\
> conda config --add channels conda-forge\
> conda install mdanalysis\
> conda install -c conda-forge rmsd\
> conda install Biopython\
> conda install argparse\
> pip install PeptideBuilder\
> conda deactivate

Make a virtual environment for Insane and martinize
> conda create -n CG_env Python=2.7\
> conda activate CG_env\
> conda install numpy\
> pip install insane\
> conda deactivate


### Requirements
Otherwise simply make sure to have installed the following modules:
- numpy
- argparse
- os
- re
- MDAnalysis
- PeptideBuilder
- Biopython
- rmsd 


and possibly insane as well:
> pip install insane

Insane is however, only compatible with python 2.7, as is the martinize script. 


# Supported python versions
Python scripts for constructing MSPs:
- Build_MSP_monomer.py - Python 3 compatible
- Assemble_two_monomers.py - Python 2+3 compatible
- Calc_number_lipids.py - Python 2+3 compatible


The Insane and Martinize scripts are however only compatible in python2.7. 
The Martinize script can also be found at the github page: https://github.com/cgmartini/martinize.py, where a python 3 version is avaliable upon request. 

### Help - documentation
For all the python scripts, additional flags and documentation can be viewed upon running the script with the flag -h.\
For the 'wrapper' scripts the flag --help can be used.

# Usage

## Constructing the atomistic MSP dimer 
Frist go into the generated virtual environment.
> conda activate ND_Builder

Use the Build_MSP_monomer.py to construct the MSP monomer given the fasta sequence.\
The output structure may appear to have clashed at the kinks between the tandem repeats. This is however solved upon a minimization.

> ./Scripts/Build_MSP_monomer.py -f base.fasta -o base

For assembling the MSP dimer, use the Assemble_two_monomers.py script

> ./Scripts/Assemble_two_monomers.py -m base.pdb -f base.fasta -w

For using default setting use the -w flag. This includes assembling the nanodisc with the LL5/5 registry, given the helix5 sequence is found in the input MSP sequence.\
A distance of 12 Å between the MSP monomers is the default setting.

For changing the registry the following flag can be used:\
--if : for selecting either LL (left-left) or RR (right-right) registry\
-t   : Changes the distance between the two MSP monomers i Angstrom. Default is 12 Å.

For changing the tandem or sequence to superimpose on the following flags can be used:\
--r1\
--r2\
-H\
When using either the --r1, --r2, or -H flags, the -f flad is needed as well, providing the fasta file.

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


## Coarse grain and insert lipids and possible a membrane protein
The script martinize and insane can then following be used for coarse graining the MSP dimer and inserting lipids and possible a membrane protein.

For insane an estimated radius of the ND is needed.\
This can be calculated with the python script

> ./Scripts/Calc_number_lipids.py -f base.fasta 

The script output a Suggestions.txt file with an estimated radius of the disc.

### Mixed lipid compositions and embedding of membrane protein 
Insane can be used for mixed lipid compositions in the ND, both with and without a possible embedded membrane protein. 

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

# Wrapper scripts for automization 
The COMMANDS_default and COMMANDS_Uni are both 'wrapper' bash scripts, which construct the nanodisc dimer using the python scripts and add lipids and solvate using insane.

The COMMANDS_default uses the deafult settings:\
Interface as Left/Left 5/5, &\
The distance between monomers as 12 Å

Within the two 'wrapper' scripts settings are set in top.\
Both scripts are for using the Martini and Charmm36 FF with the GROMACS simulation software.

> ./Wrapper_scripts/COMMANDS_Uni --help\
> ./Wrapper_scripts/COMMANDS_default --help

For using the wrapper scripts with more complex lipid compositions used the COMMANDS_Uni script and edit the option:\
lipid_flag='DLPC' in the script to eg. lipid_flag='-l POPC:90 -l POPG:10 -u POPC:90 -l POPG:10', along with the lipid_types=(POPC POPG).\
This line will be given directly to insane, which will then insert 90% POPC and 10% POPG in both the lower (-l) and upper (-u) leaflet.\
Hence the two leaflets can be constructed differently as well, by using the -l and -u flags provided to the insane script.

# Examples

## 1) Lipid only ND with simple lipid composition

Constructing 1D1 with POPC lipids.\
Go into the directory called Example1_simple_lipid_composition\
All of the below commands are embedded in the COMMANDS_default script as well.

> conda activate ND_Builder

> ./Scripts/Build_MSP_monomer.py -f 1D1.fasta -o 1D1_monomer

Assemble using default settings:

>./Scripts/Assemble_two_monomers.py -m 1D1_monomer.pdb -f 1D1.fasta -o 1D1_dimer -w

Assemble using specific settings, using tandem H5, with an RR interface and 10 angstrom space between the two MSP monomers.

>./Scripts/Assemble_two_monomers.py -m 1D1_monomer.pdb -f 1D1.fasta -o 1D1_dimer -H 'H5' -t 10 --if RR

Once the MSP dimer is constructed, lipids can be inserted with Insane, after coarse graining the MSP dimer with the martinize script. 

> conda activate CG_env

An ss.dat file is constructed, and is simply the secondary structure of the MSP defined as pure helix.\
Can be constructed automaticly as such:
> N=\`grep -c CA 1D1_dimer_assembled.pdb\`

> if test -f $"ss.dat"; then rm ss.dat; fi

> for i in \`seq 1 $N\`\
> do\
> 	echo -n 'H' >> ss.dat\
> done

> ./Scripts/martinize -f 1D1_dimer_assembled.pdb -o 1D1_cg.top -x 1D1_cg.pdb -v -name 1D1 -ss $(cat ss.dat) -ff martini22 -p Backbone

NB. If martinize comes with an error, "can't divide with zero", run the pdb through GROMACS and try again:
> gmx editconf -f 1D1_dimer_assembled.pdb -o 1D1_dimer_assembled.pdb

Next the lipids are inserted with the Insane script:\
First an estimated radius of the disc is needed.\
This can be obtained with the script:
> ./Scripts/Calc_number_lipids.py -f 1D1.fasta\
> r=\`grep 'Radius' Suggestions.txt | awk -F ':' '{print $2}'\`

> insane -f 1D1_cg.pdb -o 1D1_cg-solvent.gro -p 1D1_cg-solvent.top -pbc cubic -x 15 -y 15 -z 15 -center -l POPC -disc ${r} -a 0.7 -ring -sol W -salt 0.15 -excl -1  

> conda deactivate

Next is minimization and equlibration.\
The Run.sh or Run_ver18.sh can be used here, along with the suggested mdp files.\
This is designed for the GROMACS software.

### Circularization 

For constructing a circularized disc, the script Corr_itp_circularized.sh can be used to correct the itp file for the coarse grained MSP.\
The itp file is one of the outputs from the Martinize script.
> ./Scripts/Corr_itp_circularized.sh 1D1_A.itp

### Automization

The COMMANDS_default and COMMANDS_Uni script are wrapper scripts that can be adjusted to automize the above described process.\
For defualt settings use the COMMANDS_default script. For more specific settings, the COMMANDS_Uni script can be used.\
Settings are set within the scripts in the top.

## 2) Lipid only ND with complex lipid composition

The MSP dimer is constructed as described in the above section 'An empty ND with simple lipid composition'.\
Go into the directory called Example2_complex_lipid_composition\

> conda activate ND_Builder

> ./Scripts/Build_MSP_monomer.py -f 1D1.fasta -o 1D1_monomer

> ./Scripts/Assemble_two_monomers.py -m 1D1_monomer.pdb -f 1D1.fasta -o 1D1_dimer -w

or

>./Assemble_two_monomers.py -m 1D1_monomer.pdb -f 1D1.fasta -o 1D1_dimer -H H5 --if LL -t 12 

> N=\`grep -c CA 1D1_dimer_assembled.pdb\`\
> for i in \`seq 1 $N\`\
> do\
>       echo -n 'H' >> ss.dat\
> done

> conda activate CG_env

> ./Scripts/martinize -f 1D1_dimer_assembled.pdb -o 1D1_cg.top -x 1D1_cg.pdb -v -name 1D1 -ss $(cat ss.dat) -ff martini22 -p Backbone

> ./Scripts/Calc_number_lipids.py -f 1D1.fasta\
> r=\`grep 'Radius' Suggestions.txt | awk -F ':' '{print $2}'\`


Thereafter the insane script can easily be adjusted and used for complex lipid compositions.\
Eg. for constructing a ND with 1:9 POPG:POPC lipids, the Insane command can be like so:

> insane -f 1D1_cg.pdb -o 1D1_cg-solvent.gro -p 1D1_cg-solvent.top -pbc cubic -x 15 -y 15 -z 15 -center -l POPC:90 -l POPG:10 -u POPC:90 -l POPG:10 -disc ${r} -a 0.7 -ring -sol W -salt 0.15 -excl -1

> conda deactivate 

The -l and -u flags are for lower and upper leaflet, respectively.\
For automization, the COMMANDS_default and COMMANDS_Uni can be edited in the beginning of the file, to account for several lipid types in the disc.\
Simply change the lipid_types and lipid_flag in the top of the scripts.


## 3) A ND with a GPCR embedded - complex lipid composition 

Go into the directory called Example3_protein_embedded.\
All the commands are embedded in the COMMANDS_Uni scripts as well.


For embedding a membrane protein (MP) into a ND additional steps are required.
1) The MP is prepared, aligned, and coarse grained.
2) The MSP dimer is constructed with the python scripts, as instructed above.
3) Both the coarse grained MP and MSP are each centered in the desired box using eg. the GROMACS software.
4) The MP is inserted into the MSP dimer, by simply combining the two centered pdb files, created above.
5) Insane can then be used on the MP_MSP complex to insert the lipids, solvate, and neutralize.

For automization purposes, the COMMANDS_default and COMMANDS_Uni scripts contain the possiblity to set the MP input file and name, for embedding in the ND.\
The GCGR_cg.pdb file is avaliable here on the side for testing purposes. 
