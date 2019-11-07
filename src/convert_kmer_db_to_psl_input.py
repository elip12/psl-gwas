from large_file_processor import main, write_list, parse_args, load_pickle

# converts lines of the form 'KMER s1, s2, s3...' to lines of the form
# 's1 KMER'
# 's2 KMER'
# ...
# and writes to new file
def process(data, sim, kim, outfile):
    kmer_sample_chunk = []
    for line in data:
        linelist = line.split()
        kmer = linelist[0]
        kmer_sample_lines = []
        for sample in linelist[1:]:
            kmer_sample_lines.append(f'{sim[sample]}\t{kim[kmer]}\t1.0')
        kmer_sample_chunk.append('\n'.join(kmer_sample_lines))
    write_list(kmer_sample_chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 16
    INPUT_FILE = 'data/preprocess/kmer_sample_map_reduced.txt'
    sim_file = 'data/preprocess/sample_int_map.pickle'
    sim = load_pickle(sim_file)
    kim_file = 'data/preprocess/kmer_int_map.pickle'
    kim = load_pickle(kim_file)
    outfile = 'data/psl/contains_sample_kmer.txt'
    main(process, NUM_WORKERS, INPUT_FILE, sim=sim, kim=kim, outfile=outfile)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

