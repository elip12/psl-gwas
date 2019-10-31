from large_file_processor import LargeFileProcessor
from multiprocessing import Queue
import numpy as np
import pandas as pd

class SimilarSampleMapCreator(LargeFileProcessor):
    def __init__(self, num_workers, fname):
        super.__init__(num_workers, fname)

    def process(data, sim, outfile, n, q): 
        sample_matrix = np.zeros((n, n))
        for line in data:
            linelist = line.split()
            for i, s1 in enumerate(linelist[1: -1]):
                for s2 in linelist[i + 1:]
                    sample_matrix[sim[s1]] = sim[s2]
                    sample_matrix[sim[s2]] = sim[s1]
        q.put(sample_matrix)

    def main_wrapper(self):
        sim_file = 'data/preprocess/sample_int_map.pickle'
        sim = self.load_pickle(sim_file)
        outfile = 'data/psl/similar_sample_sample.txt'
        n = len(sim) / 2
        q = Queue()
        self.main(sim=sim, outfile=outfile, n=num_samples, q=q)
        
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
        df.to_csv('data/psl/similar_sample_sample.txt', sep='\t',
            index=False, header=False)
        

    if __name__ == '__main__':
        NUM_WORKERS = 16
        INPUT_FILE = 'data/preprocess/kmer_sample_map_reduced.txt'
        ssmc = SimilarSampleMapCreator(NUM_WORKERS, INPUT_FILE)
        ssmc.main_wrapper()

