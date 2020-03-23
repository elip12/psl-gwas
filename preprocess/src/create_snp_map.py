from large_file_processor import main, write_list, parse_args, load_pickle, dump_pickle
from multiprocessing import Queue, Manager

def process(data, sim, outfile, n, q): 
    for line in data:
    q.put((num_kmers, sample_matrix))

def main_wrapper():
    # num processes in the pool
    NUM_WORKERS = 20

    # input data
    INPUT_FILE = 'data/intermediate/unitig_sample_map.txt'
    sim_file = 'data/intermediate/sample_int_map.pickle'
    sim = load_pickle(sim_file)
    outfile = 'data/preprocessed/snp_map.txt'
    
    # number of samples
    n = int(len(sim) / 2)
    
    # multiprocessing instances for transferring data to the main thread
    m = Manager()
    q = m.Queue()

    # chunkify INPUT_FILE into NUM_WORKERS parts, create MP Pool,
    # and have each thread run process() on the chunk, with kwargs
    main(process, NUM_WORKERS, INPUT_FILE,
        sim=sim, outfile=outfile, n=n, q=q)
    
    # combine similarity matrices from each thread into a single matrix,
    # and get total number of kmers sampled
    num_kmers = 0
    sample_matrix = np.zeros((n, n))
    while not q.empty():
        q_nkmers, q_matrix = q.get()
        num_kmers += q_nkmers
        sample_matrix += q_matrix
    np.fill_diagonal(sample_matrix, np.nan) 
    
    # scale similarity counts to values in [0, 1]
    sample_matrix /= num_kmers

    df = pd.DataFrame(sample_matrix)

    # dump to tsv file for ease of restoring, and because tsv file of similarities
    # is a common input to other mGWAS programs
    df.to_csv('data/intermediate/similarities.tsv', sep='\t')
    
    # optionally read csv for ease of restoring
    #df = pd.read_csv('data/intermediate/similarities.tsv', sep='\t', header=False)
    
    #print(df)
    
    # create similarity histogram and save it
    plt.hist(df.values, bins=10, facecolor='green')
    plt.savefig('data/intermediate/hist.png', dpi=150)
    
    # set threshold
    thresh = 0.2

    # cut off all values below (less similar than) threshold
    df = df[df < thresh]
    # rescale data to [0,1]
    df *= 1/thresh
    
    # change format from
    #    c1 c2 
    # r1 a  b
    # r1 c  d
    #
    # to
    # r1 c1 a
    # r1 c2 b
    # r2 c1 c
    # r2 c2 d
    # as that is the format of PSL data files
    df = df.stack()
    df = df.reset_index()

    # write to csv 
    df.to_csv(outfile, sep='\t', index=False, header=False)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

