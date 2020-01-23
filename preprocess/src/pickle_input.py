# Reads in raw fastq data from contigs, assembles it into a dict with filenames as keys
# and lists of contigs as values, and pickles that dict.

# Usage: `python3 pickle_input.py <directory containing FASTQ files>`

import os
import pickle
from sys import argv

PICKLE_FILE = 'data/intermediate/raw.pickle'

# FASTQ files have metadata on odd lines and sequences on even lines
# ignore all errors in fastq file structure
def pickle_input(in_dir):
    seqs = {}
    for dirname, _, filenames in os.walk(in_dir):
        for filename in filenames:
            base, ext = os.path.splitext(filename)
            if ext == '.fa':
                with open(os.path.join(dirname, filename)) as f:
                    seqs[base] = [line.rstrip() for line in f \
                        if line[0] in ['A','T','C','G']]
    
    with open(PICKLE_FILE, 'wb') as pickle_file:
        pickle.dump(seqs, pickle_file)


if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError(
            'Usage: python3 pickle_input.py <directory containing FASTQ files>')
    in_dir = str(argv[1])
    if not os.path.isdir(in_dir):
        raise ValueError(
            'Usage: python3 pickle_input.py <directory containing FASTQ files>')    

    pickle_input(in_dir)

