from utility import process_file, write_dict, parse_args, printd, \
check_outfile, get_params
from multiprocessing import Queue, Manager
from collections import Counter
from preprocess_helpers import get_num_lines, create_unitig_sample_map, \
parse_input

def create_disp_nodisp_dfs(phenos):
    df = pd.read_csv(phenos, delimiter='\t')
    idcol = df.columns[0]
    df.set_index(idcol, inplace=True)
    df = df.apply(lambda x: 1.0 if x > 1.0)
    dfdisp = df.fillna(0)
    dfnodisp = 1 - df.fillna(1)
    return dfdisp, dfnodisp

def main():
    params = get_params()
    project = params['project']
    # num processes in the pool
    NUM_WORKERS = params['threads']
    # length of kmer
    k = params['k']

    # name of samples file
    samples = params['samples']
    samples = f'{project}/data/raw/{samples}'
    
    # name of phenos file
    phenos = params['phenos']
    phenos = f'{project}/data/raw/{phenos}'
    
    # threshold for num samples for each kmer, pheno pair to keep
    thresh = params['thresh']

    # dfs holding samples that display vs not display pheno
    dfdisp, dfnodisp = create_disp_nodisp_dfs(phenos)

    # read in all sequences in input into python object
    seqs = parse_input(samples)

    # number of samples
    num_samples = num_samples(samples)

    # upper and lower bounds for frequency of samples to filter kmers by
    upper = int(params['upperfreq'] * num_samples)
    lower = int(params['lowerfreq'] * num_samples)

    # input data: all input fasta files concatendated in any order
    # output data: a file containing all kmers in the population and their counts
    INPUT_FILE = f'{project}/data/preprocessed/unique_kmers.txt'
    # multiprocessing queue for transferring data to the main thread
    q = Manager().Queue()
    # chunkify INPUT_FILE into NUM_WORKERS parts, create MP Pool,
    # and have each thread run process() on the chunk, with kwargs
    process_file(create_unitig_sample_map, NUM_WORKERS, INPUT_FILE,
        raw=seqs, q=q, k=k, thresh=thresh, upper=upper, lower=lower,
        dfdisp=dfdisp, dfnodisp=dfnodisp)
    
    outfile = f'{project}/data/preprocessed/unitig_sample_map.txt'
    check_outfile(outfile)
    # write all chunks to output file sequentially
    while not q.empty():
        write_list(q.get(), outfile)
   

if __name__ == '__main__':
    parse_args()
    main()

