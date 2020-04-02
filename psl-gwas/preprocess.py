#!/usr/bin/env python3
from utility import process_file, write_dict, parse_args, printd, \
file_exists, get_params
from multiprocessing import Queue, Manager
from collections import Counter
from preprocess_helpers import num_samples, create_unitig_sample_map, \
parse_input, similar_sample
import int_maps
from os.path import join

# separate phenos df into 2 dfs, one holding phenos present and 1 holding
# phenos absent.
def create_disp_nodisp_dfs(phenos):
    df = pd.read_csv(phenos, delimiter='\t')
    idcol = df.columns[0]
    df.set_index(idcol, inplace=True)
    df = df.apply(lambda x: 1.0 if x > 1.0)
    df = df.apply(lambda x: 0.0 if x < 0.0)
    dfdisp = df.fillna(0)
    dfnodisp = 1 - df.fillna(1)
    return dfdisp, dfnodisp

def main():
    # get params
    params = get_params()
    project = params['project']
    k = params['k']
    samples = params['samples']
    phenos = params['phenos']
    thresh = params['thresh']

    # define file paths
    unique_kmers_file = join(project, 'data', 'preprocessed', 'unique_kmers.txt')
    phenos_file = join(project, 'data', 'raw', phenos)
    samples_file = join(project, 'data', 'raw', samples)
    similarities_tsv = join(project, 'data', 'preprocessed', 'sample_similarities.tsv')
    hist_orig_file = join(project, 'data', 'preprocessed', 'hist_orig.png')
    hist_scaled_file = join(project, 'data', 'preprocessed', 'hist_scaled.png')
    similar_sample_file = join(project, 'data', 'preprocessed', 'similar_sample_sample.txt')
    unitigs_file = join(project, 'data', 'preprocessed', 'unitig_sample_map.txt')
    sim_file = join(project, 'data', 'preprocessed', 'sample_int_map.txt') 
    pim_file = join(project, 'data', 'preprocessed', 'pheno_int_map.txt')
    uim_file = join(project, 'data', 'preprocessed', 'unitig_int_map.txt')

    # create and load sample and pheno int maps
    if not file_exists(sim_file):
        int_maps.create_sim_file(samples, sim_file)
    if not file_exists(pim_file):
        int_maps.create_pim_file(phenos, pim_file)
    sim = load_pickle(sim_file)
    
    if not file_exists(unitigs_file) or not file_exists(similar_sample_file):
        # dfs holding samples that display vs not display pheno
        dfdisp, dfnodisp = create_disp_nodisp_dfs(phenos)
        # read in all sequences in input into python object
        seqs = parse_input(samples)
        # number of samples
        n_samples = num_samples(samples)
        # upper and lower bounds for frequency of samples to filter kmers by
        upper = int(params['upperfreq'] * n_samples)
        lower = int(params['lowerfreq'] * n_samples)
        # multiprocessing queue for transferring data to the main thread
        q = Manager().Queue()
        process_file(create_unitig_sample_map, unitigs_file,
            raw=seqs, q=q, k=k, thresh=thresh, upper=upper, lower=lower,
            dfdisp=dfdisp, dfnodisp=dfnodisp, sim=sim, n=n_samples)
       
        sample_matrix = np.zeros((n, n), dtype=np.uint32)
        num_kmers = 0
        # write all chunks to output file sequentially
        create_unitigs_file = False
        if not file_exists(unitigs_file):
            create_unitigs_file = True
        while not q.empty():
            unitigs, q_num_kmers, q_sample_matrix = q.get()
            if create_unitigs_file:
                write_list(unitigs, unitigs_file)
            num_kmers += q_num_kmers
            sample_matrix += q_sample_matrix
        
        if not file_exists(similar_sample_file):
            similar_sample(sample_matrix, num_kmers, similarities_tsv,
                hist_orig_file, hist_scaled_file, similar_sample_file)
    
    if not file_exists(uim_file):
        int_maps.create_unitig_int_map(unitigs_file, uim_file)

if __name__ == '__main__':
    parse_args()
    main()

