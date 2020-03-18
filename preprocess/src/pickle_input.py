# Reads in raw fastq data from contigs, assembles it into a dict with filenames as keys
# and lists of contigs as values, and pickles that dict.

import os
import pickle
from sys import argv

PICKLE_FILE = 'data/intermediate/raw.pickle'

# FASTQ files have metadata on odd lines and sequences on even lines
# ignore all errors in fastq file structure
def pickle_input(infile):
    seqs = {}
    with open(infile, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        sample, fname = line.split('\t')
        with open(fname, 'r') as f:
            contiglines = f.readlines()
        seqs[sample] = [contigline.rstrip() for contigline in contiglines \
            if contigline[0] in ['A','T','C','G']]
    
    with open(PICKLE_FILE, 'wb') as pickle_file:
        pickle.dump(seqs, pickle_file)

if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError(
            'Usage: python3 pickle_input.py <samples>.tsv')
    infile = argv[1]
    pickle_input(in_dir)

