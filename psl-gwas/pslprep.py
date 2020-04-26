#!/usr/bin/env python3
###############################################################################
##  pslprep.py
##  This script creates the files PSL takes as input from other preprocessed
##  and raw data.
###############################################################################
from utility import (process_file, write_list, parse_args, load_pickle, printd,
    get_params, file_exists)
from pslprep_model import (unitig_sample_db, unitig_pheno_db, sample_pheno,
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
    uim_file = join(project, 'data', 'preprocessed', 'unitig_int_map.pkl')
    unitig_sample_map_file = join(project, 'data', 'preprocessed', 'unitig_sample_map.txt')
    unitig_pheno_map_file = join(project, 'data', 'preprocessed', 'unitig_pheno_map.txt')
    phenos_file = join(project, 'data', 'raw', params['phenos'])
    contains_sample_unitig_file = join(project, 'data', 'preprocessed', 'contains_sample_unitig.txt')
    value_sample_pheno_file = join(project, 'data', 'preprocessed', 'value_sample_pheno.txt')
    value_unitig_pheno_file = join(project, 'data', 'preprocessed', 'value_unitig_pheno.txt')
    similar_pheno_pheno_file = join(project, 'data', 'preprocessed', 'similar_pheno_pheno.txt')
    
    sim = load_pickle(sim_file)
    pim = load_pickle(pim_file)
    
    # simulate truth data
    if params.get('truth'):
        truths_infile = join(project, 'data', 'raw', params['truth'])
        truths_dict = create_truths_dict(truths_infile, pim)
        truth_unitig_pheno_file = join(project, 'data', 'preprocessed', 'truth_unitig_pheno.txt')
    else:
        truths_dict = None
        truth_unitig_pheno_file = None

    # create smaller psl input files that can be efficiently done w 1 thread
    if not file_exists(value_sample_pheno_file):
        sample_pheno(phenos_file, sim, pim, value_sample_pheno_file)
    if not file_exists(similar_pheno_pheno_file):
        similar_pheno(phenos_file, pim, similar_pheno_pheno_file)
    
    contains_exists = file_exists(contains_sample_unitig_file)
    value_exists = file_exists(value_unitig_pheno_file)
    truths_exists = file_exists(truth_unitig_pheno_file) if params.get('truth') else True
    
    lock = Manager().Lock()
    
    if not contains_exists:
        process_file(unitig_sample_db, unitig_sample_map_file,
            uim_file=uim_file, lock=lock, truths=truths_dict,
            contains_sample_unitig_file=contains_sample_unitig_file)

    if not value_exists or not truths_exists:
        process_file(unitig_pheno_db, unitig_pheno_map_file,
            uim_file=uim_file, value_unitig_pheno_file=value_unitig_pheno_file,
            truth_unitig_pheno_file=truth_unitig_pheno_file,
            lock=lock, truths=truths_dict) 

if __name__ == '__main__':
    parse_args()
    main()

