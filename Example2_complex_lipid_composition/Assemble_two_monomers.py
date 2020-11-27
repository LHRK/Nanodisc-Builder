#!/usr/bin/env python

import numpy as np
import MDAnalysis as md
import argparse
import os

parser = argparse.ArgumentParser(description='Script for assembling two monomers of a MSP for a Nanodisc. The flag -w can be set for default settings. Otherwise distance between monomers (-t), the tandem sequence for superimposing (-H) or just the two residues for superimposing the two monomers into a dimer (--res1, --res2) can be set. If the default settings with the flas -w does not work, you might have mutations in the tandem sequence reffered to as Helix5.')

parser.add_argument('-m', dest='monomer', action='store', type=str, help='-m for monomers pdb file (Required)')      
parser.add_argument('-f', dest='fasta', action='store', type=str, help='-f for fasta file (Required)')      
parser.add_argument('-o', dest='final', action='store', type=str, help='-o for finale output pdb name (Required)')
parser.add_argument('-w', action='store_true', help='set for default settings (Optional)')

parser.add_argument('-H', dest='Helix', action='store', type=str, help='-H for which tandem repeats to use for superimposing the dimers. Default=H5 (Optional)', default=None)
parser.add_argument('-t', dest='trans', action='store', type=int, help='-t for the space between the two monomers in angstrom (Optional)', default=None)
parser.add_argument('--if', dest='interface', action='store', type=str, help='--if for interface, either LL or RR (Stands for Left-Left and Right-Right) (Optional)', default=None)
parser.add_argument('--r1', dest='res1', action='store', type=int, help='-r1 for residue or starting residue to superimpose with (Optional)', default=None)
parser.add_argument('--r2', dest='res2', action='store', type=int, help='-r2 for end residue to superimpose with (Optional)', default=None)                                                                                  
                                                                                                       
args = parser.parse_args()                                                                             

def check_interface(u, interface='LL'):
    sele  = u.select_atoms('name CA')
    L = []
    R = []
    for idx, res in enumerate(sele.resnames):
        if res == 'PRO':
            left = [idx+1, idx+4, idx+8]
            right = [idx+3, idx+6, idx+10]
            L.append(left)
            R.append(right)
    
    Left_resids = " ".join(["{0:d}".format(i) for j in L for i in j])
    Right_resids = " ".join(["{0:d}".format(i) for j in R for i in j])
    L_coor = u.select_atoms('resid {0:s}'.format(Left_resids)).positions
    R_coor = u.select_atoms('resid {0:s}'.format(Right_resids)).positions
    
    H = np.mean(u.select_atoms('all').positions, axis=0)
    
    L_c = H - np.mean(L_coor, axis=0)
    R_c = H- np.mean(R_coor, axis=0)
    L_c_nor = L_c / np.linalg.norm(L_c)    
    R_c_nor = R_c / np.linalg.norm(R_c)    
    L_dot = np.dot(L_c_nor, np.array([0,0,1]))
    R_dot = np.dot(R_c_nor, np.array([0,0,1]))
    if interface == 'LL':
        if L_dot > 0:
            return_value = 1
        else:
            return_value = -1
    elif interface=='RR':
        if R_dot > 0:
            return_value = 1
        else:
            return_value = -1
    return return_value

#################################################################################################################
## Tandem sequences are from:                                                                                  ##
## Grinkova, Yelena V., Ilia G. Denisov, and Stephen G. Sligar.                                                ##
## "Engineering extended membrane scaffold proteins for self-assembly of soluble nanoscale lipid bilayers."    ##
## Protein Engineering, Design and Selection 23.11 (2010): 843-848.                                            ##
#################################################################################################################

def auto_identify (fasta, H='H5'):
	Tandems    = {"H1":"LKLLDNWDSVTSTFSKLREQLG",
	              "H2":"PVTQEFWDNLEKETEGLRQEMS",
	              "H3":"KDLEEVKAKVQ",
	              "H4":"PYLDDFQKKWQEEMELYRQKVE",
	              "H5":"PLRAELQEGARQKLHELQEKLS",
	              "H6":"PLGEEMRDRARAHVDALRTHLA",
	              "H7":"PYSDELRQRLAARLEALKENGG",
	              "H8":"ARLAEYHAKATEHLSTLSEKAK",
	              "H9":"PALEDLRQGLL",
	              "H10":"PVLESFKVSFLSALEEYTKKLNTQ",
	}
	
	Helix = Tandems[H]
	fasta_in = open(fasta).readlines()[1:]
	fasta_seq = " ".join([i.strip('\n') for i in fasta_in])
	begin_seq = fasta_seq.find(Helix)
	end_seq = begin_seq + len(Helix)
	return begin_seq, end_seq


if args.w == True:
	res1, res2 = auto_identify (args.fasta, H='H5')
	interface='LL'
	trans = 12
else:
	res1 = args.res1
	res2 = args.res1
	interface = args.interface
	trans = args.trans
	Helix_sel = args.Helix

if args.Helix != None:
	res1, res2 = auto_identify (args.fasta, H=args.Helix)

if args.trans != None:
	trans = args.trans

u = md.Universe(args.monomer)
prot = u.select_atoms('all').positions
helix = u.select_atoms('resid {0:d}-{1:d}'.format(res1, res2))
H = helix.positions

C = np.mean(H, axis= 0 )
cosa, sina, = C[:2] / np.sqrt(C[0]**2 + C[1]**2)
Rz = np.array([ [cosa, -sina, 0], [sina, cosa, 0], [0, 0, 1] ])
Hz = np.dot(prot-np.mean(prot, axis=0), Rz)
u.select_atoms('all').positions = Hz
out1 = u.select_atoms('all')
out1.write('{0:s}'.format(args.final+'_monomer1.pdb'))

a = np.radians(180)
#Ry = np.array([ [np.cos(a), 0, np.sin(a)], [0, 1, 0], [-np.sin(a), 0, np.cos(a)] ])
Rx = np.array([ [1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)] ])

inter_corr = check_interface(u, interface)

Z = np.array([ [0,0,trans*inter_corr] ])

#Hzy = np.dot(Hz-np.mean(Hz, axis=0)+Z, Ry.T)
Hzx = np.dot(Hz-np.mean(Hz, axis=0)+Z, Rx)

# +Z or -Z depending on the interface

#u.select_atoms('all').positions = Hzy
u.select_atoms('all').positions = Hzx
out = u.select_atoms('all')
out.write('{0:s}'.format(args.final+'_monomer2.pdb'))

os.system("sed -i 's/ S / A /g' {0:s}".format(args.final+'_monomer1.pdb'))
os.system("sed -i 's/ A / B /g' {0:s}".format(args.final+'_monomer2.pdb'))
os.system("cat {0:s} {1:s} > {2:s}".format(args.final+'_monomer1.pdb', args.final+'_monomer2.pdb', args.final+'_assembled.pdb'))

os.system("sed -i '/END/d' {0:s}".format(args.final+'_assembled.pdb'))
