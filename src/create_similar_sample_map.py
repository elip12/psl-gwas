from large_file_processor import main, write_list, parse_args, load_pickle
from multiprocessing import Queue, Manager
import numpy as np
import pandas as pd
import random

def process(data, sim, outfile, n, q): 
    sample_matrix = np.zeros((n, n))
    for line in data:
        if random.random() > 0.1:
            continue
        linelist = line.split()
        for i, s1 in enumerate(linelist[1: -1]):
            for s2 in linelist[i + 1:]:
                sample_matrix[int(sim[s1])][int(sim[s2])] += 1 
                sample_matrix[int(sim[s2])][int(sim[s1])] += 1
    q.put(sample_matrix)

def main_wrapper():
    NUM_WORKERS = 32
    INPUT_FILE = 'data/preprocess/kmer_sample_map_reduced.txt'
    sim_file = 'data/preprocess/sample_int_map.pickle'
    sim = load_pickle(sim_file)
    outfile = 'data/psl/similar_sample_sample.txt'
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
    df = pd.DataFrame(sample_matrix)
    df = df.stack()
    df = df.reset_index()
    df.to_csv(outfile, sep='\t', index=False, header=False)

if __name__ == '__main__':
    parse_args()
    main_wrapper()
