#!/usr/bin/env python

import numpy as np
import argparse


parser = argparse.ArgumentParser(description='Scirpt for determine a radius of the ND and the number of lipids given a certain area per lipid in squared angstrom')

parser.add_argument('-f', dest='fasta', action='store', type=str, help='-f for fasta file input')

args = parser.parse_args()
fasta = args.fasta 

#simply deleting the first line in fasta file
f = open('{0:s}'.format(fasta), 'r').readlines()
f_out = open('seq.txt', 'w')
for line in f[1:]:
    f_out.write('{0:s}'.format(line))
f_out.close()

file_seq = ''.join(open('seq.txt').read().split())


nres          = len(file_seq)
C             = (nres / 3.6) * 1.5
r             = (nres * 1.5) / (2*np.pi)


file_suggest = open('Suggestions.txt', 'w')
file_suggest.write('Radius (nm) of the ND  :{0:.2f}\n'.format((r/10)))
file_suggest.close()
 
