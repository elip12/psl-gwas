from utility import process_file, write_dict, parse_args, printd, \
check_outfile, get_params
from multiprocessing import Queue, Manager
from collections import Counter
from preprocess import reduce_input

def main():
    params = get_params()
    # num processes in the pool
    NUM_WORKERS = params['threads']

    # length of kmer
    k = params['k']

    # input data: all input fasta files concatendated in any order
    # output data: a file containing all kmers in the population and their counts
    INPUT_FILE = 'data/intermediate/unique_kmers.txt'
    check_outfile(outfile)

    # multiprocessing queue for transferring data to the main thread
    q = Manager().Queue()

    # chunkify INPUT_FILE into NUM_WORKERS parts, create MP Pool,
    # and have each thread run process() on the chunk, with kwargs
    process_file(process, NUM_WORKERS, INPUT_FILE,
        q=q, k=k)
    
    counter = Counter()
    while not q.empty():
        counter.update(q.get())
    printd('Finished consolidating counters.')
   
    # write counter to file
    write_dict(counter, outfile)

if __name__ == '__main__':
    parse_args()
    main()

