from large_file_processor import LargeFileProcessor

class KmerDbCleaner(LargeFileProcessor):
    def __init__(self, num_workers, fname):
        super.__init__(num_workers, fname)

    def process(data, kmer_db_file, THRESH):
        kmer_db_chunk = []
        for line in data:
            linelist = line.split()
            if len(linelist) < THRESH + 1:
                continue
            kmer_db_line = [linelist[0]]
            for csv in linelist[1:]:
                values = csv.split(',')
                kmer_db_line.append(values[0])
            kmer_db_chunk.append(' '.join(line_to_write))
        self.write_list(kmer_db_chunk, kmer_db_file)

    def main_wrapper(self):
        kdf = 'data/preprocess/kmer_sample_map_reduced.txt'
        THRESH = 10
        self.main(kmer_db_file=kdf, THRESH=THRESH)

    if __name__ == '__main__':
        NUM_WORKERS = 16
        INPUT_FILE = 'data/preprocess/kmer_sample_map.txt'
        kdc = KmerDbCleaner(NUM_WORKERS, INPUT_FILE)
        kdc.main_wrapper()

