#!/usr/bin/env python3
###############################################################################
##  count_kmers.py
##  Extracts all unique kmers from input data, and maps them to the number
##  of times they occur.
###############################################################################
from utility import process_file, write_dict, parse_args, printd, \
file_exists, get_params, write_list, complement
from multiprocessing import Queue, Manager
from collections import Counter
from os import remove
from os.path import join
from int_maps import create_unitig_int_map

# cat all sample .fa together
def cat_samples(samples, outfile):
    with open(samples, 'r') as f:
        lines = f.readlines() 
    for line in lines:
        sample, fname = line.rstrip().split('\t')
        if not (fname.endswith('.fa') or fname.endswith('.fsa')):
            continue
        with open(fname, 'r') as f:
            write_list(f.readlines(), outfile)

# takes in a chunk of contigs, and creates a counter holding all the kmers
# and their counts in that chunk.
def process(data, q, k):
    counter = Counter()
    testset = set('ATCG')
    for line in data:
        if line.startswith('>') or len(line) < k:
            continue
        for i in range(len(line) - k):
            kmer = line[i: i + k]
            if set(kmer) != testset: # kmers with Ns and Ys in them
                continue
            comp = complement(kmer)
            counter[min(kmer, comp)] += 1
    q.put(counter)

def main():
    # load params 
    params = get_params()
    project = params['project']
    k = params['k']

    # define file paths
    samples_file = join(project, 'data', 'raw', params['samples'])
    outfile = join(project, 'data', 'preprocessed', 'unique_kmers.txt')
    catted_samples = join(project, 'data', 'preprocessed', 'samples.fa')

    # check if output file exists; if so, do nothing.
    if file_exists(outfile):
        exit(0)

    # create catted samples file if it does not exist.
    if not file_exists(catted_samples):
        cat_samples(samples_file, catted_samples)

    # multiprocessing queue for transferring data to the main thread
    q = Manager().Queue()

    # invoke process(...) on catted_samples files with kwargs, for each thread
    process_file(process, catted_samples, q=q, k=k)
    
    # consolidate all threads' counters into single counter holding all kmers
    counter = Counter()
    while not q.empty():
        counter.update(q.get())
    printd('Finished consolidating counters.')

    # write counter to file
    write_dict(counter, outfile, sep='\t')
    
    # remove catted samples file
    if file_exists(catted_samples):
        remove(catted_samples)

if __name__ == '__main__':
    parse_args()
    main()

