#!/usr/bin/env python3
###############################################################################
##  preprocess.py
##  This script maps each kmer to the samples it appears in, takes a random
##  sample of kmers and constructs a sample similarity matrix, consolidates
##  kmers into longer unitigs, filters unitigs for correlation with pheno,
##  and creates integer mappings for samples, phenos, and unitigs to
## speed up the association test.
###############################################################################
from utility import process_file, parse_args, printd, \
file_exists, get_params, load_pickle, write_list
from multiprocessing import Queue, Manager
from collections import Counter
from preprocess_model import num_samples, create_unitig_sample_map, \
parse_input, similar_sample, create_disp_nodisp_dfs
import int_maps
from os.path import join
import numpy as np

def main():
    # get params
    params = get_params()
    project = params['project']

    # define file paths
    unique_kmers_file = join(project, 'data', 'preprocessed', 'unique_kmers.txt')
    phenos_file = join(project, 'data', 'raw', params['phenos'])
    samples_file = join(project, 'data', 'raw', params['samples'])
    similarities_tsv = join(project, 'data', 'preprocessed', 'sample_similarities.tsv')
    hist_orig_file = join(project, 'data', 'preprocessed', 'hist_orig.png')
    hist_scaled_file = join(project, 'data', 'preprocessed', 'hist_scaled.png')
    similar_sample_file = join(project, 'data', 'preprocessed', 'similar_sample_sample.txt')
    unitigs_file = join(project, 'data', 'preprocessed', 'unitig_sample_map.txt')
    sim_file = join(project, 'data', 'preprocessed', 'sample_int_map.pkl') 
    pim_file = join(project, 'data', 'preprocessed', 'pheno_int_map.pkl')
    uim_file = join(project, 'data', 'preprocessed', 'unitig_int_map.pkl')

    # create and load sample and pheno int maps
    if not file_exists(sim_file):
        int_maps.create_sample_int_map(samples_file, phenos_file, sim_file)
    if not file_exists(pim_file):
        int_maps.create_pheno_int_map(phenos_file, pim_file)
    sim = load_pickle(sim_file)
    
    # only do processing if output files do not exist
    if not file_exists(unitigs_file) or not file_exists(similar_sample_file):
        # dfs holding samples that display vs not display pheno
        dfdisp, dfnodisp = create_disp_nodisp_dfs(phenos_file, sim)
        # read in all sequences in input into python object
        seqs = parse_input(samples_file)
        # number of samples
        n_samples = num_samples(samples_file)
        # upper and lower bounds for frequency of samples to filter kmers by
        upper = int(0.95 * n_samples)
        lower = int(0.01 * n_samples)
        # multiprocessing queue for transferring data to the main thread
        q = Manager().Queue()
        process_file(create_unitig_sample_map, unique_kmers_file,
            raw=seqs, q=q, k=params['k'], thresh=params['thresh'], upper=upper,
            lower=lower, dfdisp=dfdisp, dfnodisp=dfnodisp, sim=sim, n=n_samples)
       
        sample_matrix = np.zeros((n_samples, n_samples))
        num_kmers = 0
        # write all chunks to output files sequentially
        create_unitigs_file = False
        if not file_exists(unitigs_file):
            create_unitigs_file = True
        while not q.empty():
            unitigs, q_num_kmers, q_sample_matrix = q.get()
            if create_unitigs_file:
                write_list(unitigs, unitigs_file)
            num_kmers += q_num_kmers
            sample_matrix += q_sample_matrix
        
        # create sample similarity file
        if not file_exists(similar_sample_file):
            similar_sample(sample_matrix, num_kmers, similarities_tsv,
                hist_orig_file, hist_scaled_file, similar_sample_file)
    # create unitig int map
    if not file_exists(uim_file):
        int_maps.create_unitig_int_map(unitigs_file, uim_file)

if __name__ == '__main__':
    parse_args()
    main()

