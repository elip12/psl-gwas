from large_file_processor import LargeFileProcessor

class RKCMCreator(LargeFileProcessor):
    def __init__(self, num_workers, fname):
        super.__init__(num_workesr, fname)

    def process(data, cim, rkc_file):
        rkc_chunk = []
        for line in data:
            linelist = line.split()
            rkc_line = [f'{linelist[0]} {cim[c]}' for c in cim]
            rkc_chunk.append('\n'.join(rkc_line))
        self.write_list(rkc_chunk, rkc_file)

    def main_wrapper(self):
        cim_file = 'data/preprocess/class_int_map.pickle'
        cim = self.load_pickle(cim_file)
        rkc = 'data/psl/resistance_kmer_class.txt'
        self.main(cim=cim, rkc_file=rkc)

    if __name__ == '__main__':
        NUM_WORKERS = 16
        INPUT_FILE = 'data/preprocess/kmer_sample_map_reduced.txt'
        rkcm = RKCMCreator(NUM_WORKERS, INPUT_FILE)
        rkcm.main_wrapper()

