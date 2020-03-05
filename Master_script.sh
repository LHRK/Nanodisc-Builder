#!/bin/bash

## COMMENTS ###
# I had to add some antifreez particles to the system


# You will only need this script and the fasta file for the Nanodisc you want to build. 
# The script will copy the nescarry files to the working directory
# It build the MSPs atomistically with the python scripts, and then coarse grain the structures and insert lipids using insane.py

# Stuff you need to edit
export base='2n5e' #Base name for the system / Nanodisc
FF='martini22' #martini22 or martini_3.0.b.3.2/martini303v.partition
fasta="${base}.fasta" 
r1=144 #Resids in the sequence for imposing
r2=144
trans_z=12 #distance in angstrom between the two MSP monomers
interface='RR' #LL = Left-Left interface or RR for Right-Right interface. 
x_len=15 #Size of the box in nm
y_len=15
z_len=15
lipid_type='DLPC' # lipid type. The script insane insert the lipids
r=3.3 # Radius of the phosphorlipid bilayer in the nanodisc
apl=0.59 # Area per lipid in nm² for use in insane
asym_corr=-3 # Use for correting the amount of lipids in each leaflet, if insane dose not inseret the same amount in each leaflet. 
script_dir='/home/au447022/Documents/GU/ND/Building/SCIPTS' # Directory for where all the scripts, force field parameters and mdp files are located. 

########################
# Reads in flags
Standard=''
var=$1

f=${var:-Standard}


if [ "$f" = "Standard" ]
then
echo 'Sligar disc, with terminals'
        elif [ "$f" = "Connect" ]
        then
        echo 'Circularized, terminals connected'
	cp ${script_dir}/Corr_itp_circularized.sh .
                elif [ "$f" = "--help" ]
                then
                echo -ne "Edit input parameters in this scirpt and run using the flag Connect for circularized, or none for Standard discs with terminals.\nFor cicularized discs, the itp file needs modification.\nThe correction of the itp file for the circularized is currently only supported for Coarse Grained.\n"
		exit
                        elif [ "$f" != "Connect" ]
                        then
                        echo 'Use Connect for circularized or none for standard sligar discs'
fi
########################
#Setting the water name depending on the martini verison

if [ "$FF" = "martini22" ]
then
W_name="W"
cp  ${script_dir}/martini_3.0.b.3.2/martinize .
cp ${script_dir}/Martini2_itpFiles/*.itp .
	elif [ "$FF" != "martini22" ]
	then
	W_name="WN"
        cp -r ${script_dir}/martini_3.0.b.3.2 .
        cp ${script_dir}//martini_3.0.b.3.2/martinize .
fi

########################

# Building the MSPs and assembling them atomistically
cp ${script_dir}/Build_MSP_monomer.py .
cp ${script_dir}/Assemble_two_monomers.py .

source activate py2.7
python -W ignore Build_MSP_monomer.py -f $fasta -o ${base}_monomer
python -W ignore Assemble_two_monomers_updated.py -m ${base}_monomer.pdb -o ${base}_dimer --r1 ${r1} --r2 ${r2} -t ${trans_z} --if ${interface}
source deactivate

#Adding hydrogens and generating an atomistic topology - so far only for sligar, so for circularized you will need to correct the itp file
gmx editconf -f ${base}_dimer_assembled.pdb -o ${base}-box.gro -box ${x_len} ${y_len} ${z_len}
cp -r ${script_dir}/charmm36.ff .
gmx pdb2gmx -f ${base}-box.gro -o ${base}-pdb2gmx.gro -p ${base}_AA_pdb2gmx.top -ff charmm36 -water none -i ${base}_posres.itp


line=`grep -n 'endif' ${base}_AA_pdb2gmx.top | gawk -F : '{print $1}'`
head -n $((line)) ${base}_AA_pdb2gmx.top > ${base}_AA.itp
sed -i '20d' ${base}_AA.itp
sed -i '20d' ${base}_AA.itp


tot_lines=`wc -l ${base}_AA_pdb2gmx.top | gawk '{print $1}'`
top_lines=$((tot_lines-line))
tail -n $top_lines ${base}_AA_pdb2gmx.top > ${base}_AA.top


sed -i "s/Protein/${base}/" ${base}_AA.top
sed -i "s/Protein/${base}/" ${base}_AA.itp
sed -i "1i #include \"charmm36.ff/forcefield.itp\"" ${base}_AA.top
sed -i "2i #include \"${base}_AA.itp\"" ${base}_AA.top

cp ${base}_AA.top ${base}_AA_nowat.top


sed -i "3i #include \"$lipid_type.itp\"" ${base}_AA.top
sed -i "4i #include \"charmm36.ff/tip3p.itp\"" ${base}_AA.top
sed -i "5i #include \"charmm36.ff/ions.itp\"" ${base}_AA.top

gmx editconf -f ${base}_dimer_assembled.pdb -box ${x_len} ${y_len} ${x_len} -o ${base}_dimer_box.pdb
AA_struc="${base}_dimer_box.pdb"

#Generating the coarse grained structure
source activate py2.7
##Generate the ss.dat file
N=`grep -c CA ${AA_struc}`
rm ss.dat # Just in case there is one file already
for i in `seq 1 $N`
do
echo -n 'H' >> ss.dat
done

################ Coarse graining ######################################
if [ "$f" = "Connect" ]
then
	echo "./martinize -f ${AA_struc} -o ${base}_cg.top -x ${base}_cg.pdb -v -name ${base} -ss $(cat ss.dat) -ff ${FF} -p Backbone -nt"
	./martinize -f ${AA_struc} -o ${base}_cg.top -x ${base}_cg.pdb -v -name ${base} -ss $(cat ss.dat) -ff ${FF} -p Backbone -nt
elif [ "$f" = "Standard" ]
then
	echo "./martinize -f ${AA_struc} -o ${base}_cg.top -x ${base}_cg.pdb -v -name ${base} -ss $(cat ss.dat) -ff ${FF} -p Backbone"
	./martinize -f ${AA_struc} -o ${base}_cg.top -x ${base}_cg.pdb -v -name ${base} -ss $(cat ss.dat) -ff ${FF} -p Backbone
fi
source deactivate

#######################################################################

#Correction of the itp file for circularized discs, currently only for coarse grained itp files
if [ "$f" = "Connect" ]
then
	./Corr_itp_circularized.sh ${base}_A.itp
fi

############# Adding lipids, water and NaCl ##########################
source activate custom
echo "insane -f ${base}_cg.pdb -o ${base}_cg-solvent.gro -p ${base}_cg-solvent.top -pbc cubic -x ${x_len} -y ${y_len} -z ${z_len} -center -l ${lipid_type} -disc ${r} -a ${apl} -ring -sol ${W_name} -salt 0.15 -excl -1 -asym ${asym_corr}"
insane -f ${base}_cg.pdb -o ${base}_cg-solvent.gro -p ${base}_cg-solvent.top -pbc cubic -x ${x_len} -y ${y_len} -z ${z_len} -center -l ${lipid_type} -disc ${r} -a ${apl} -ring -sol ${W_name} -salt 0.15 -excl -1 -asym ${asym_corr}
source deactivate

#####################################################################

# For adding it to the AA top file
count_lipids=0
for i in `grep "$lipid_type " ${base}_cg-solvent.top | gawk '{print $2}'`
do
count_lipids=$((count_lipids+i))
done

line=`grep -n "$base " ${base}_AA.top | gawk -F ':' '{print $1}'`
line=$((line+1))

sed -i "${line}i $lipid_type     $count_lipids" ${base}_AA.top

#Correcting the names of the ions in the gro file
#updating the top file with the correct names and itp files

if [ "$FF" = "martini22" ]
then 
sed -i 's/"martini.itp"/"martini_v2.2.itp"/' ${base}_cg-solvent.top
sed -i "2i #include \"${base}_A.itp\"" ${base}_cg-solvent.top
sed -i "4i #include \"${lipid_type}.itp\"" ${base}_cg-solvent.top
sed -i '5i #include "martini_v2.0_solvents.itp"' ${base}_cg-solvent.top
sed -i '6i #include "martini_v2.0_ions.itp"' ${base}_cg-solvent.top
sed -i "s/Protein/${base}_A/" ${base}_cg-solvent.top
sed -i "15i ${base}_A          1" ${base}_cg-solvent.top

	elif [ "$FF" != "martini22" ]
	then
	sed -i 's/NA+/TNA/g' ${base}_cg-solvent.gro
	sed -i 's/CL-/TCL/g' ${base}_cg-solvent.gro
	sed -i 's/"martini.itp"/"martini_3.0.b.3.2\/martini_v3.0.b.3.2.itp"/' ${base}_cg-solvent.top
	sed -i "2i #include \"${base}_A.itp\"" ${base}_cg-solvent.top
	sed -i '4i #include "martini_3.0.b.3.2/martini_v3.0_phospholipids.itp"' ${base}_cg-solvent.top
	sed -i '5i #include "martini_3.0.b.3.2/martini_v3.0_solvents.itp"' ${base}_cg-solvent.top
	sed -i '6i #include "martini_3.0.b.3.2/martini_v3.0_ions.itp"' ${base}_cg-solvent.top
	sed -i "s/Protein/${base}_A/" ${base}_cg-solvent.top
	sed -i "15i ${base}_A          1" ${base}_cg-solvent.top
	sed -i 's/NA+/TNA/' ${base}_cg-solvent.top
	sed -i 's/CL-/TCL/' ${base}_cg-solvent.top
fi

cp ${script_dir}/MDP_FILES_CG/*.mdp .

#Grompping a min tpr
gmx grompp -f min.mdp -c ${base}_cg-solvent.gro -o ${base}_cg-solvent-min.tpr -p ${base}_cg-solvent.top -maxwarn 1

##Generating a index file with the group solvent being water + ions
echo 'q' | gmx make_ndx -f ${base}_cg-solvent-min.tpr -o ${base}_cg-solvent.ndx 
sel1=`grep '\[' ${base}_cg-solvent.ndx | gawk '{print NR-1, $2}' | grep "${W_name}" | gawk '{print $1}'`
sel2=`grep '\[' ${base}_cg-solvent.ndx | gawk '{print NR-1, $2}' | grep 'ION' | gawk '{print $1}'`
echo "$sel1 | $sel2" > sel
echo "q" >> sel
cat sel | gmx make_ndx -n ${base}_cg-solvent.ndx -o ${base}_cg-solvent.ndx

sel4=`grep '\[' ${base}_cg-solvent.ndx | gawk '{print NR-1, $2}' | grep "${W_name}_ION" | gawk '{print $1}'`
echo "name $sel4 solvent" > sel5
echo "q" >> sel5

cat sel5 | gmx make_ndx -n ${base}_cg-solvent.ndx -o ${base}_cg-solvent.ndx
rm sel sel5

cp ${script_dir}/Run.sh .

############# Minimazing and equlibrating CG ##########################
#./Run.sh
#rm step*.pdb
#rm \#*#


