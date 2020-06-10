from utility import process_file, write_list, parse_args, load_pickle, \
get_params, file_exists, write_files
from multiprocessing import Manager
from os.path import join

def process(data, lock, pim, kim_file, thresh, fsa_file, scored_kmers_file):
    kim = load_pickle(kim_file)
    chunk = []
    for line in data:
        linelist = line.split()
        if float(linelist[2]) < thresh:
            continue
        outline = (kim[int(linelist[0])], pim[int(linelist[1])], linelist[2])
        chunk.append(outline)
    kmers = [f'>{i}\n{line[0]}' for i, line in enumerate(chunk)]
    values = ['\t'.join(tup) for tup in chunk]
    write_files(lock,
        (values, scored_kmers_file),
        (kmers, fsa_file))

def main():
    # get params
    params = get_params()
    project = params['project']
    
    # define file paths
    INPUT_FILE = join(project, 'data', 'postprocessed', 'KMERPHENO.txt')
    pim_file = join(project, 'data', 'preprocessed', 'pheno_int_map.pkl')
    fsa_file = join(project, 'data', 'postprocessed', 'scored_kmers.fsa')
    kim_file = join(project, 'data', 'preprocessed', 'kmer_int_map.pkl')
    scored_kmers_file = join(project, 'data', 'postprocessed', 'scored_kmers.txt')

    # create output files if they do not exist
    if file_exists(fsa_file):
        fsa_file = None
    if file_exists(scored_kmers_file):
        scored_kmers_file = None
    if fsa_file or scored_kmers_file:
        lock = Manager().Lock()
        pim = load_pickle(pim_file)
        
        process_file(process, INPUT_FILE, lock=lock, pim=pim, kim_file=kim_file,
                thresh=params['classification-thresh'], fsa_file=fsa_file,
                scored_kmers_file=scored_kmers_file)

if __name__ == '__main__':
    parse_args()
    main()

