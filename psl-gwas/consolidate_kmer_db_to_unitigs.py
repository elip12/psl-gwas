from large_file_processor import main, write_list, parse_args

def process(data, outfile, k):
    prev_line = data[0]
    prev_linelist = prev_line.split()
    prev_unitig = prev_linelist[0]
    prev_samples = prev_linelist[1:]
    prev_samples = [','.join(s.split(',')[:2]) for s in prev_samples]
    unitig_db_chunk = []
    for line in data[1:]:
        linelist = line.split()
        this_unitig = linelist[0]
        this_samples = linelist[1:]
        this_samples = [','.join(s.split(',')[:2]) for s in this_samples]
        
        # if kmers are sequential and the same set of samples contain both kmers
        if prev_unitig[-(k - 1):] == line[0:k - 1] \
                and len(this_samples) == len(prev_samples) \
                and set(this_samples) == set(prev_samples):
            this_unitig = f'{prev_unitig}{line[k - 1]}'
            prev_line = f'{this_unitig}{prev_line[len(prev_unitig):]}'
        else:
            unitig_db_chunk.append(prev_line)
            prev_line = line
        
        prev_unitig = this_unitig
        prev_samples = this_samples
    write_list(unitig_db_chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 20
    K = 30
    INPUT_FILE = 'data/intermediate/kmer_sample_map.txt'
    outfile = 'data/intermediate/unitig_sample_map.txt'
    main(process, NUM_WORKERS, INPUT_FILE, outfile=outfile, k=K)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

