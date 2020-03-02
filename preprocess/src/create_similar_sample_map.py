from large_file_processor import main, write_list, parse_args, load_pickle, dump_pickle
from multiprocessing import Queue, Manager
import numpy as np
import pandas as pd
import random
from matplotlib import pyplot as plt


def process(data, sim, outfile, n, q): 
    sample_matrix = np.zeros((n, n))
    for line in data:
        if random.random() > 0.001:
            continue
        linelist = line.split()
        for i, s1_ in enumerate(linelist[1: -1]):
            s1 = s1_.split(',')[0]
            for s2_ in linelist[i + 1:]:
                s2 = s2_.split(',')[0]
                sample_matrix[int(sim[s1])][int(sim[s2])] += 1 
                sample_matrix[int(sim[s2])][int(sim[s1])] += 1
    q.put(sample_matrix)

def main_wrapper():
    NUM_WORKERS = 16
    INPUT_FILE = 'data/intermediate/kmer_sample_map.txt'
    sim_file = 'data/intermediate/sample_int_map.pickle'
    sim = load_pickle(sim_file)
    outfile = 'data/preprocessed/similar_sample_sample.txt'
    n = int(len(sim) / 2)
    m = Manager()
    q = m.Queue()
    main(process, NUM_WORKERS, INPUT_FILE,
        sim=sim, outfile=outfile, n=n, q=q)
    
    sample_matrix = np.zeros((n, n))
    while not q.empty():
        sample_matrix += q.get()
    np.fill_diagonal(sample_matrix, np.nan)
    a_min = np.nanmin(sample_matrix)
    a_max = np.nanmax(sample_matrix)
    
    # scale to value in [0, 1]
    sample_matrix = (sample_matrix - a_min) / (a_max - a_min)
    dump_pickle(sample_matrix, 'data/intermediate/relatedness_matrix.pickle')
    df = pd.DataFrame(sample_matrix)

    df.to_csv('similarities.tsv', sep='\t')
    print(df)
    plt.hist(df.values, bins=10)
    plt.savefig('hist.png', dpi=150)
    '''
        create histogram
        drop everything below threshhold
    '''
    thresh = 0.2
    df = df[df < thresh]
    
    df = df.stack()
    df = df.reset_index()
    df.to_csv(outfile, sep='\t', index=False, header=False)

if __name__ == '__main__':
    parse_args()
    main_wrapper()
