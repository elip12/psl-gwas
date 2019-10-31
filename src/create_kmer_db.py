from large_file_processor import LargeFileProcessor

class KmerDbCreator(LargeFileProcessor):
    def __init__(self, num_workers, fname):
        super.__init__(num_workers, fname)

    # DNA base complement
    def complement(self, base):
        if base == 'A':
            return 'T'
        elif base == 'T':
            return 'A'
        elif base == 'G':
            return 'C'
        else:
            return 'G'

    # Creates a dictionary of all kmers passed in data, and their complements
    # Then, iterates through genomes. If a genome contains a kmer, that genome's
    # metadata is added to dict entry for that kmer.
    # After all samples have been iterated over, writes kmers to file.
    def process(self, data, raw, K, outfile):
        kmers = {}
        for line in data:
            kmer = line.split(' ')[0]
            comp = ''.join(map(complement, kmer.split()))
            kmers[kmer] = ''
            kmers[comp] = ''
        
        for count, (raw_id, seq) in enumerate(raw.items()):
            for c_id, contig in enumerate(seq):
                l = len(contig)
                if l >= K: # ensure this contig is long enough to sample
                    for i in range(l - K + 1):
                        kmer = contig[i: i + K]
                        if kmer in kmers:
                            kmers[kmer] += f' {raw_id},{c_id},{i}'
            self.printd(f'\tProcessed genome {count + 1}', end='\r')
        self.write_dict(kmers, outfile)
    
    # defines some variables we need to pass to process, then calls main
    def main_wrapper():
        pickle_file = 'data/preprocess/raw.pickle'
        raw = self.load_pickle(pickle_file)
        K = 30
        outfile = 'data/preprocess/kmer_sample_map.txt'
        self.main(raw=raw, K=k, outfile=outfile)

if __name__ == '__main__':
    NUM_WORKERS = 16
    INPUT_FILE = 'data/preprocess/unique_kmers_reduced.txt'
    kdc = KmerDbCreator(NUM_WORKERS, INPUT_FILE)
    kdc.main_wrapper()
