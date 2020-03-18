from large_file_processor import main, write_list, parse_args, load_pickle

def process(data, cim, outfile):
    kim_file = 'data/intermediate/kmer_int_map.pickle'
    kim = load_pickle(kim_file)
    chunk = []
    for line in data:
        linelist = line.split()
        if float(linelist[2]) < 0.35: #0.9999999:
            continue
        outline = f'{kim[int(linelist[0])]}\t{cim[int(linelist[1])]}\t{linelist[2]}'
        chunk.append(outline)
    write_list(chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 20
    INPUT_FILE = 'inferred-predicates/KMERRESISTANCE.txt'
    cim_file = 'data/intermediate/class_int_map.pickle'
    cim = load_pickle(cim_file)
    outfile = 'data/postprocessed/scored_kmers.txt'
    main(process, NUM_WORKERS, INPUT_FILE, cim=cim, outfile=outfile)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

