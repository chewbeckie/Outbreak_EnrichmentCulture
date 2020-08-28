#!/usr/bin/env python3
import sys

gfa = open(sys.argv[1], 'r')
fasta = open(sys.argv[2], 'w')

for line in gfa:
	if line[0] != 'S': continue
	dat = line.rstrip().split('\t')
	fasta.write('>'+dat[1]+'\n'+dat[2]+'\n')

