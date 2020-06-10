#!/usr/bin/env python3
###############################################################################
##  pslprep.py
##  This script creates the files PSL takes as input from other preprocessed
##  and raw data.
###############################################################################
from utility import (process_file, write_list, parse_args, load_pickle, printd,
    get_params, file_exists)
from pslprep_model import (kmer_sample_db, kmer_pheno_db, sample_pheno,
    similar_pheno, create_truths_dict)
from multiprocessing import Manager
from os.path import join

def main():
    # get params
    params = get_params()
    project = params['project']

    # define data paths
    sim_file = join(project, 'data', 'preprocessed', 'sample_int_map.pkl')
    pim_file = join(project, 'data', 'preprocessed', 'pheno_int_map.pkl')
    kim_file = join(project, 'data', 'preprocessed', 'kmer_int_map.pkl')
    kmer_sample_map_file = join(project, 'data', 'preprocessed', 'kmer_sample_map.txt')
    kmer_pheno_map_file = join(project, 'data', 'preprocessed', 'kmer_pheno_map.txt')
    phenos_file = join(project, 'data', 'raw', params['pheno'])
    contains_sample_kmer_file = join(project, 'data', 'preprocessed', 'contains_obs.txt')
    value_sample_pheno_file = join(project, 'data', 'preprocessed', 'samplePheno_obs.txt')
    value_kmer_pheno_file = join(project, 'data', 'preprocessed', 'kmerPheno_target.txt')
    similar_pheno_pheno_file = join(project, 'data', 'preprocessed', 'similarPheno_obs.txt')
    
    sim = load_pickle(sim_file)
    pim = load_pickle(pim_file)
    
    # incorporate truth data
    if params.get('truth'):
        truths_infile = join(project, 'data', 'raw', params['truth'])
        truths_dict = create_truths_dict(truths_infile, pim)
        truth_kmer_pheno_file = join(project, 'data', 'preprocessed', 'kmerPheno_truth.txt')
    else:
        truths_dict = None
        truth_kmer_pheno_file = None

    # incorporate baseline data
    if params.get('baseline'):
        baseline_infile = join(project, 'data', 'raw', params['baseline'])
        baseline_dict = create_truths_dict(baseline_infile, pim)
        baseline_kmer_pheno_file = join(project, 'data', 'preprocessed', 'baseline_obs.txt')
    else:
        baseline_dict = None
        baseline_kmer_pheno_file = None

    # create smaller psl input files that can be efficiently done w 1 thread
    if not file_exists(value_sample_pheno_file):
        sample_pheno(phenos_file, sim, pim, value_sample_pheno_file)
    if not file_exists(similar_pheno_pheno_file):
        similar_pheno(phenos_file, pim, similar_pheno_pheno_file)
    
    contains_exists = file_exists(contains_sample_kmer_file)
    value_exists = file_exists(value_kmer_pheno_file)
    truths_exists = file_exists(truth_kmer_pheno_file) if params.get('truth') else True
    baseline_exists = file_exists(baseline_kmer_pheno_file) if params.get('baseline') else True

    lock = Manager().Lock()
    
    if not contains_exists:
        process_file(kmer_sample_db, kmer_sample_map_file,
            kim_file=kim_file, lock=lock, truths=truths_dict,
            contains_sample_kmer_file=contains_sample_kmer_file)

    if not value_exists or not truths_exists or not baseline_exists:
        if value_exists:
            value_kmer_pheno_file = None
        if truths_exists:
            truth_kmer_pheno_file = None
        if baseline_exists:
            baseline_kmer_pheno_file = None
        process_file(kmer_pheno_db, kmer_pheno_map_file,
            kim_file=kim_file, value_kmer_pheno_file=value_kmer_pheno_file,
            truth_kmer_pheno_file=truth_kmer_pheno_file,
            lock=lock, truths=truths_dict, baseline=baseline_dict,
            baseline_kmer_pheno_file=baseline_kmer_pheno_file) 

if __name__ == '__main__':
    parse_args()
    main()

