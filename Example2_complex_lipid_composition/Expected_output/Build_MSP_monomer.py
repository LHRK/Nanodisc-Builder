#!/usr/bin/env python

import numpy as np
import PeptideBuilder
from PeptideBuilder import Geometry
import MDAnalysis as mda
import re
import argparse
import rmsd
import Bio.PDB
import os

parser = argparse.ArgumentParser(description='Scirpt for building a monomer of a MSP given a fasta sequence.')

parser.add_argument('-f', dest='fasta', action='store', type=str, help='-f for fasta file input')
parser.add_argument('-o', dest='final', action='store', type=str, help='-o for finale output pdb name')

args = parser.parse_args()

def seq_extract (input_seq, output='One_letter_seq.txt'):
    """Is reads in a fasta file and writes out the 3-letter
    sequence to an ourput file and returns the 3-letter sequence."""
    file = open(input_seq)
    lines = file.readlines()
    list = []
    for l in lines:
            if re.match('^>', l):
                    continue
            else:
                    for a in l:
                            if re.match('^\n', a):
                                    continue
                            else:
                                    list.append(a)
    seq_out = open(output+'_three_letter_seq.txt', 'w')
    seq = []
    mapping1_3 = {"A":"ALA",
              "R":"ARG",
              "N":"ASN",
              "D":"ASP",
              "C":"CYS",
              "Q":"GLN",
              "E":"GLU",
              "G":"GLY",
              "H":"HIS",
              "I":"ILE",
              "L":"LEU",
              "K":"LYS",
              "F":"PHE",
              "M":"MET",
              "P":"PRO",
              "S":"SER",
              "T":"THR",
              "W":"TRP",
              "Y":"TYR",
              "V":"VAL"}
    for a in list:
        seq_out.write('{0:s}\n'.format(mapping1_3[a]))
        seq.append(mapping1_3[a])
    seq_out.close()
    return seq



def CA_trace_coor(seq, omega):
    """The function generates a CA trace in a circle
    depending on the input sequence lenght.
    It returns the coordinates"""
    rad = 5.5 # helix radius, selected from article mentioned below. Use this radius for calculate the area of the phospholipid patch in the disc
    t   = 100 # helix turn angle
    l   = 1.5 # translation   #### 1.6 Angstrom in translation fit well with the homology model built
    
    nres = len(seq)
    coor = np.zeros([nres,3])
    R0 = int(np.ceil((nres * l) / (2 * np.pi)))
    
    ########################################################################################################################
    #  The equation for the radius (R0) and the values for the helix radius r is from:                                     #
    #  Denisov, I. G. et al. "Directed self-assembly of monodisperse phosphorlipid bilayer Nanodiscs with controlled size" #
    #  Journal of the American Chemical Society 126.11 (2004): 3477-3487                                                   #
    #  The translation l for an ideal alpha helix is 1.5 Angstrom                                                          #
    ########################################################################################################################
    
    nturns = int(np.ceil(nres / 4)) + 1 # I have added an additional turn as a buffer, to make sure the first and last residues does not overlap
    dt_alpha = 2 * np.pi /nturns
    t = np.radians(t)
    
    res = 1
    
    for it in range(nturns+1):
        alpha = it * dt_alpha
        Orx = R0 * np.cos(alpha)
        Ory = R0 * np.sin(alpha)
        Orz = 0
        for  r in range(4):
            a = r * t + omega
            x = Orx + r * np.sin(a) * np.cos(alpha) - l * r * np.sin(alpha)
            y = Ory + r * np.sin(a) * np.sin(alpha) + l * r * np.cos(alpha)
            z = Orz - r * np.cos(a)
            try:
                coor[res-2,0] = x
                coor[res-2,1] = y
                coor[res-2,2] = z
            except IndexError:
                 #print ('Last residue was {0:d}'.format(res))
                 break
            res = res + 1
    return coor


def build_entire_helix(file_seq):
    """This function builds one straight
    alpha-helix from the input sequence"""
    geo = Geometry.geometry(file_seq[0])
    geo.phi=-60
    geo.psi_im1 =-40
    structure = PeptideBuilder.initialize_res(geo)
    for r in file_seq[1:]:
        structure = PeptideBuilder.add_residue(structure , r, -60, -40)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save( 'Entire_helix.pdb')
    return

def align_kabsch(a1, a2, helix, trace):
    """This function aligns each helix piece with a corresponding CA-trace.
    The striaght helix is divided into helix pieces where either a Proline residue is present
    or where there is double Glycine (GG)"""
    Q = trace
    Q_gec = np.mean(Q, axis=0)
    H = Helix.select_atoms('resid {0:d}-{1:d}'.format(a1,a2)).positions
    gec = np.mean(H, axis = 0)
    P = Helix.select_atoms('resid {0:d}-{1:d} and name CA'.format(a1,a2)).positions - gec # translate
    Q = Q-Q_gec
    R = rmsd.kabsch(P, Q) #calculate the rotations matrix
    H2 = H - gec # Translate the entire helix
    coors = np.dot(H2, R) + Q_gec # Rotate and translate back
    return coors 


def rotate_apolar(a1, a2, name):
    """The function translate and rotate the helix piece,
    such that the apolar side of the helix is towards
    the centrum of the disc"""
    Apolar = 'ALA VAL ILE LEU MET PHE TYR TRP PRO'
    #Apolar = 'ALA VAL ILE LEU MET PHE TYR TRP'
    u = mda.Universe(name)
    helix = u.select_atoms('resid {0:d}-{1:d}'.format(a1,a2))
    H = helix.positions
    C = np.mean(helix.positions, axis=0) #Centrum point of the helix
    apolar_coor = u.select_atoms('resname {0:s} and resid {1:d}-{2:d}'.format(Apolar, a1, a2)).positions
    apolar = np.mean(apolar_coor, axis=0)
    P = C - apolar
    cosa, sina = C[:2]/np.sqrt(C[0]**2 + C[1]**2) #Only considering x and y values
    #Rotation matrix around z
    Rz = np.array([ [cosa, -sina, 0], [sina, cosa, 0], [0, 0, 1] ])
    Hz = np.dot(H - C, Rz) # Translate and rotate around z the entire helix 
    Pz = np.dot(P, Rz) #Rotate the hydrophobic vector around z as well
    cosb, dummy, sinb = Pz / np.sqrt(Pz[0]**2 + Pz[2]**2)
    Ry = np.array([ [cosb, 0, -sinb], [0, 1, 0], [sinb, 0, cosb] ])
    Hy = (np.dot( np.dot( Hz, Ry), Rz.T) + C)
    #Here we rotate around y to get the invers hydrophobic vector P along x axis
    #And then we 'back' rotate around z and undo the translation
    #helix.positions = Hy
    #helix.write(out)
    return Hy

fasta = args.fasta 
seq = seq_extract(fasta)


#simply deleting the first line in fasta file
f = open('{0:s}'.format(fasta), 'r').readlines()
f_out = open('seq.txt', 'w')
for line in f[1:]:
    f_out.write('{0:s}'.format(line))
f_out.close()

file_seq = ''.join(open('seq.txt').read().split())

coordinats = CA_trace_coor(seq, 1)
build_entire_helix(file_seq)

Helix = mda.Universe('Entire_helix.pdb')

count = 0
C = []
for idx, s in enumerate(file_seq):
    if s == 'P':
        C.append(align_kabsch(count+1, idx, Helix, coordinats[np.arange(count,idx), :]))
        count = idx
    elif file_seq[0] == 'G':
        continue
    elif file_seq[-1] == 'G':
        break 
    elif s =='G' and file_seq[idx+1]=='G':
        C.append(align_kabsch(count+1, idx, Helix, coordinats[np.arange(count,idx), :]))
        count = idx 
    else:
        continue
C.append(align_kabsch(count+1, len(file_seq), Helix, coordinats[np.arange(count,len(file_seq)), :]))

Helix.select_atoms('all').positions = np.vstack((C))
Helix.select_atoms('all').write('Helix_circular.pdb')


count = 0
AP = []
for idx, s in enumerate(file_seq):
    if s == 'P':
        AP.append(rotate_apolar(count+1, idx, 'Helix_circular.pdb'))
        count = idx 
    elif file_seq[0] == 'G':
        continue
    elif file_seq[-1] == 'G':
        break
    elif s =='G' and file_seq[idx+1]=='G':
        AP.append(rotate_apolar(count+1, idx, 'Helix_circular.pdb'))
        count = idx 
    else:
        continue
AP.append(rotate_apolar(count+1, len(file_seq), 'Helix_circular.pdb'))

Helix_AP = mda.Universe('Helix_circular.pdb')
Helix_AP.select_atoms('all').positions = np.vstack((AP))
Helix_AP.select_atoms('all').write(args.final+'.pdb')

os.system("rm Entire_helix.pdb One_letter_seq.txt_three_letter_seq.txt seq.txt Helix_circular.pdb")



