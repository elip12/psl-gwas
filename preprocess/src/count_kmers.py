from large_file_processor import main, write_list, parse_args, load_pickle, dump_pickle, printd
from multiprocessing import Queue, Manager
import numpy as np
import pandas as pd
import random
from matplotlib import pyplot as plt
from collections import Counter

def process(data, q, k):
    counter = Counter()
    for line in data:
        if line.startswith('>') or len(line) < k:
            continue
        for i in range(len(line) - k):
            kmer = line[i: i + k]
            counter[kmer] += 1
    q.put(counter)

def main_wrapper():
    # 4 processes in the pool
    NUM_WORKERS = 4

    # length of kmer
    k = 30

    # input data: all input fasta files concatendated in any order
    INPUT_FILE = 'data/intermediate/samples.fa'
    # output data: a file containing all kmers in the population and their counts
    outfile = 'data/intermediate/unique_kmers.txt'
    
    # multiprocessing instances for transferring data to the main thread
    m = Manager()
    q = m.Queue()

    # chunkify INPUT_FILE into NUM_WORKERS parts, create MP Pool,
    # and have each thread run process() on the chunk, with kwargs
    main(process, NUM_WORKERS, INPUT_FILE,
        q=q, k=k)
    
    counter = Counter()
    while not q.empty():
        counter += q.get()
    printd('Finished consolidating counters.')
    with open(outfile, 'w') as f:
        f.writelines(f'{k}\t{v}\n' for k, v in counter.items())
    printd('Finished writing to outfile.')

if __name__ == '__main__':
    parse_args()
    main_wrapper()

