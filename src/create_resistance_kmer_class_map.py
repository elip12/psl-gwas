from large_file_processor import main, write_list, parse_args, load_pickle

def process(data, cim, outfile):
    kim_file = 'data/preprocess/kmer_int_map.pickle'
    kim = load_pickle(kim_file)
    rkc_chunk = []
    for line in data:
        linelist = line.split()
        lc = int(len(cim) / 2)
        rkc_line = [f'{kim[linelist[0]]} {c}' for c in range(lc)]
        rkc_chunk.append('\n'.join(rkc_line))
    write_list(rkc_chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 16
    INPUT_FILE = 'data/preprocess/kmer_sample_map_reduced.txt'
    cim_file = 'data/preprocess/class_int_map.pickle'
    cim = load_pickle(cim_file)
    outfile = 'data/psl/resistance_kmer_class.txt'
    main(process, NUM_WORKERS, INPUT_FILE, cim=cim, outfile=outfile)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

