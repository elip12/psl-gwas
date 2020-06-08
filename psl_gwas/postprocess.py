from utility import process_file, write_list, parse_args, load_pickle, \
get_params, file_exists, write_files
from multiprocessing import Manager
from os.path import join

def process(data, lock, pim, uim_file, fsa_file, scored_unitigs_file):
    uim = load_pickle(uim_file)
    chunk = []
    for line in data:
        linelist = line.split()
        outline = (uim[int(linelist[0])], pim[int(linelist[1])], linelist[2])
        chunk.append(outline)
    unitigs = [f'>{i}\n{line[0]}' for i, line in enumerate(chunk)]
    values = ['\t'.join(tup) for tup in chunk]
    write_files(lock,
        (values, scored_unitigs_file),
        (unitigs, fsa_file))

def main():
    # get params
    params = get_params()
    project = params['project']
    
    # define file paths
    INPUT_FILE = join(project, 'data', 'postprocessed', 'UNITIGPHENO.txt')
    pim_file = join(project, 'data', 'preprocessed', 'pheno_int_map.pkl')
    fsa_file = join(project, 'data', 'postprocessed', 'scored_unitigs.fsa')
    uim_file = join(project, 'data', 'preprocessed', 'unitig_int_map.pkl')
    scored_unitigs_file = join(project, 'data', 'postprocessed', 'scored_unitigs.txt')

    # create output files if they do not exist
    if file_exists(fsa_file):
        fsa_file = None
    if file_exists(scored_unitigs_file):
        scored_unitigs_file = None
    if fsa_file or scored_unitigs_file:
        lock = Manager().Lock()
        pim = load_pickle(pim_file)
        
        process_file(process, INPUT_FILE, lock=lock, pim=pim, uim_file=uim_file,
                fsa_file=fsa_file, scored_unitigs_file=scored_unitigs_file)

if __name__ == '__main__':
    parse_args()
    main()

