#!/usr/bin/env python3 
from os import listdir
from os.path import basename, splitext, join

'''
annotate psl-gwas kmers output with prodigal output.
gets the contig/gene each kmer appears in by cross referencing 
'''

dirname = 'prodigal'
for fname in listdir(dirname):
    bname = basename(fname) # get file name without path
    bname, ext = splitext(bname)
    if ext != '.fna':
        continue
    with open(join(dirname, fname), 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('>'):
            newline = line.rstrip()
            newline += f' # {bname}\n'
            lines[i] = newline
    with open(join(dirname, fname), 'w') as f:
        f.writelines(lines)

