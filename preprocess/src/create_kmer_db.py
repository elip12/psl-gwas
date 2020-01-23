from large_file_processor import main, printd, write_dict, parse_args, \
    load_pickle

DEBUG = False

# DNA base complement
def complement(base):
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
def process(data, raw, K, outfile):
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
        printd(f'\tProcessed genome {count + 1}', end='\r')
    write_dict(kmers, outfile)
    
# defines some variables we need to pass to process, then calls main
def main_wrapper():
    NUM_WORKERS = 16
    INPUT_FILE = 'data/intermediate/unique_kmers_reduced.txt'
    pickle_file = 'data/intermediate/raw.pickle'
    raw = load_pickle(pickle_file)
    K = 30
    outfile = 'data/intermediate/kmer_sample_map.txt'
    main(process, NUM_WORKERS, INPUT_FILE, raw=raw, K=K, outfile=outfile)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

