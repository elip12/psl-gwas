from large_file_processor import LargeFileProcessor

class KmerDbToPSLInputConverter(LargeFileProcessor):
    def __init__(self, num_workers, fname):
        super.__init__(num_workers, fname)

    # converts lines of the form 'KMER s1, s2, s3...' to lines of the form
    # 'KMER s1'
    # 'KMER s2'
    # ...
    # and writes to new file
    def process(self, data, sim, outfile):
        kmer_sample_chunk = []
        for line in data:
            linelist = line.split()
            kmer = linelist[0]
            kmer_sample_lines = []
            for sample in linelist[1:]:
                kmer_sample_lines.append(f'{kmer}\t{sim[sample]}\t1.0')
            kmer_sample_chunk.append('\n'.join(kmer_sample_lines))
        self.write_list(kmer_sample_chunk, outfile)

    def main_wrapper(self):
        sim_file = 'data/preprocess/sample_int_map.pickle'
        sim = self.load_pickle(sim_file)
        outfile = 'data/psl/contains_sample_kmer.txt'
        self.main(sim=sim, outfile=outfile)

if __name__ == '__main__':
    NUM_WORKERS = 16
    INPUT_FILE = 'data/preprocess/sample_kmer_map_reduced.txt'
    k = KmerDbToPSLInputConverter(NUM_WORKERS, INPUT_FILE)
    k.main_wrapper()

