#!/bin/bash

### Bash script for backmapping the Nanodisc, re-solvate and neutralize using gromacs followed by a two step minimization
### If perfered, the initram.sh script can also be applied

scirpt_dir='PATH/SCRIPTS' #directory containing the scripts and mdp files

cp ${script_dir}/MDP_FILES_AA/*.mdp .
cp -r ${script_dir}/Mapping .
cp ${script_dir}/backward.py .


############### STUFF TO EDIT ###################################
                                                                #
base='2N5E' #Base name of the system                            #
lipid_type='DLPC' # in the coarse grained structure for grep    #
lipid_type_AA='DMPC' # The lipid type for backmapping           #
                                                                #
# The dimensions for the Box                                    #
x_len=15                                                        #
y_len=15                                                        #
z_len=13                                                        #
## NB                                                           #
#It is assumed that the CG input file is named base_cg.gro      #
                                                                #
#################################################################

gmx editconf -f ${base}_cg.gro -o ${base}_cg.pdb
# Grep only protein and membrane
grep -e BB -e SC -e ${lipid_type} ${base}_cg.pdb > ${base}_cg_nowat.pdb 

# Backmap CG structure to AA
./backward.py -f ${base}_cg_nowat.pdb -p ${base}_AA.top -o ${base}_AA.gro 

# Resolvate and neutralize
echo '1' | gmx editconf -f ${base}_AA.gro -box ${x_len} ${y_len} ${z_len} -o ${base}_AA_box.gro -princ yes

cp ${base}_AA.top ${base}_AA_wat.top

#Use a custom vdwradii.dat file where the radius for C and P are increased, such that water are not placed in the bilayer

gmx solvate -cp ${base}_AA_box.gro -p ${base}_AA_wat.top -o ${base}_AA_wat.gro 

gmx grompp -f min.mdp -c ${base}_AA_wat.gro -p ${base}_AA_wat.top -o ${base}_AA_wat.tpr -maxwarn 1 

cp ${base}_AA_wat.top ${base}_AA_solv.top

echo 'q' | gmx make_ndx -f ${base}_AA_wat.tpr -o ${base}_AA_wat.ndx
grep '\[' ${base}_AA_wat.ndx | gawk '{print NR-1, $2}' | grep SOL | gawk '{print $1}' > sel

cat sel | gmx genion -s ${base}_AA_wat.tpr -p ${base}_AA_solv.top -o ${base}_AA_solv.gro -conc 0.15 -neutral
rm sel

#Generation of index file
echo 'q' | gmx make_ndx -f ${base}_AA_solv.gro -o ${base}_AA_solv.ndx
grep '\[' ${base}_AA_solv.ndx | gawk '{print NR-1, $2}' | grep ${lipid_type_AA} | gawk '{print $1}' > temp
m=`head -n 1 temp`
echo "name $m Membrane" > sel
s=`grep '\[' ${base}_AA_solv.ndx | gawk '{print NR-1, $2}' | grep Water_and_ions | gawk '{print $1}'`
echo "name $s Solvent " >> sel
echo 'q' >> sel

cat sel | gmx make_ndx -f ${base}_AA_solv.gro -o ${base}_AA_solv.ndx
rm sel
rm temp

## The below mdp files are from the initram.sh, the number of steps in the two minimization files are just double up

# Step 1 mini
gmx grompp -f 1-EM.mdp -c ${base}_AA_solv.gro -n ${base}_AA_solv.ndx -p ${base}_AA_solv.top -o 1-EM.tpr
gmx mdrun -v -deffnm 1-EM

# Step 2 mini 
gmx grompp -f 2-EM.mdp -c 1-EM.gro -n ${base}_AA_solv.ndx -p ${base}_AA_solv.top -o 2-EM.tpr
gmx mdrun -v -deffnm 2-EM

# Step 3 MD_1 
gmx grompp -f 3-mdpr-0.0002.mdp -c 2-EM.gro -p ${base}_AA_solv.top -n ${base}_AA_solv.ndx -o 3-mdpr-0.0002.tpr -r 2-EM.gro
gmx mdrun -v -deffnm 3-mdpr-0.0002
