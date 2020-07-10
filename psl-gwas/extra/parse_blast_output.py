#!/usr/bin/env python3 
from sys import argv

'''
after you blast your kmers against some database, use this to get a file
of the genes/contigs each kmer is associated with
'''

kmersfile = argv[1]
infile = argv[2]
outfile = argv[3]
with open(kmersfile, 'r') as f:
    kmers = f.readlines()
with open(infile, 'r') as f:
    lines = f.readlines()
with open(outfile, 'w') as f:
    for line in lines:
        if line.startswith('Query='):
            kmernum = int(line.split('=')[1].rstrip())
            kmer = kmers[kmernum * 2 + 1]
            f.write(kmer)
        if line.startswith('>'):
            f.write(line)

