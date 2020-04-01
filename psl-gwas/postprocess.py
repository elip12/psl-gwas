from utility import process_file, write_list, parse_args, load_pickle, \
get_params
from multiprocessing import Manager, Queue

def process(data, q, pim, uim_file):
    uim = load_pickle(uim_file)
    chunk = []
    for line in data:
        linelist = line.split()
        if float(linelist[2]) < 0.95:
            continue
        outline = (uim[int(linelist[0])], pim[int(linelist[1])], linelist[2])
        chunk.append(outline)
    q.put(chunk)

def main():
    INPUT_FILE = f'{project}/data/postprocessed/KMERRESISTANCE.txt'
    pim_file = f'{project}/data/preprocessed/pheno_int_map.pkl'
    fsa_file = f'{project}/data/postprocessed/scored_unitigs.fsa'
    scored_unitigs_file = f'{project}/data/postprocessed/scored_unitigs.txt'
    
    params = get_params()
    project = params['project']
    NUM_WORKERS = params['threads']
    check_outfile(outfile)
    q = Manager().Queue()
    pim = load_pickle(pim_file)
    
    process_file(process, NUM_WORKERS, INPUT_FILE, q=q, pim=pim, uim_file=uim_file)
    
    while not q.empty():
        chunk = q.get()
        unitigs = [f'>{i}\n{line[0]}\n' for i, line in enumerate(chunk)]
        write_list(unitigs, fsa_file)
        write_list(chunk, scored_unitigs_file)

if __name__ == '__main__':
    parse_args()
    main()

