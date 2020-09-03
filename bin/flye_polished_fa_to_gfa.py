#!/usr/bin/env python3
import sys

gfa_in = open(sys.argv[1], 'r')
fasta_in = open(sys.argv[2], 'r')
gfa_out = open(sys.argv[3], 'w')

# copy in the polished sequences
for line in fasta_in:
	if(line[0]=='>'):
		gfa_out.write('S\t'+line.rstrip()[1:])
	else:
		gfa_out.write('\t'+line)

# copy the links from the other gfa
for line in gfa_in:
	if(line[0]!='L'): continue
	gfa_out.write(line)

