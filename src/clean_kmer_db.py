from large_file_processor import main, write_list, parse_args

def process(data, outfile, THRESH):
    kmer_db_chunk = []
    for line in data:
        linelist = line.split()
        if len(linelist) < THRESH + 1:
            continue
        kmer_db_line = [linelist[0]]
        for sci in linelist[1:]:
            values = sci.split(',')
            kmer_db_line.append(values[0])
        kmer_db_chunk.append(' '.join(kmer_db_line))
    write_list(kmer_db_chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 16
    INPUT_FILE = 'data/preprocess/kmer_sample_map.txt'
    outfile = 'data/preprocess/kmer_sample_map_reduced.txt'
    THRESH = 10
    main(process, NUM_WORKERS, INPUT_FILE, outfile=outfile, THRESH=THRESH)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

