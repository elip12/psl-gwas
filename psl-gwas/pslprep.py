#!/usr/bin/env python3
###############################################################################
##  pslprep.py
##  This script creates the files PSL takes as input from other preprocessed
##  and raw data.
###############################################################################
from utility import process_file, write_list, parse_args, load_pickle, printd, \
get_params, file_exists
from pslprep_model import unitig_db, sample_pheno, similar_pheno
from multiprocessing import Manager, Queue
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
    phenos_file = join(project, 'data', 'raw', params['phenos'])
    contains_sample_unitig_file = join(project, 'data', 'preprocessed', 'contains_sample_unitig.txt')
    value_sample_pheno_file = join(project, 'data', 'preprocessed', 'value_sample_pheno.txt')
    value_unitig_pheno_file = join(project, 'data', 'preprocessed', 'value_unitig_pheno.txt')
    similar_pheno_pheno_file = join(project, 'data', 'preprocessed', 'similar_pheno_pheno.txt')
    
    # simulate truth data
    if params.get('truths'):
        truths_infile = join(project, 'data', 'raw', 'genes.fa') # this might change
        truths_dict = create_truths_dict(truths_infile)
        truths_outfile = join(project, 'data', 'preprocess', 'truth_pheno_unitig.txt')
    else:
        truths_dict = None

    sim = load_pickle(sim_file)
    pim = load_pickle(pim_file)

    # create smaller psl input files that can be efficiently done w 1 thread
    if not file_exists(value_sample_pheno_file):
        sample_pheno(phenos_file, sim, pim, value_sample_pheno_file)
    if not file_exists(similar_pheno_pheno_file):
        similar_pheno(phenos_file, pim, similar_pheno_pheno_file)
    
    contains_exists = file_exists(contains_sample_unitig_file)
    value_exists = file_exists(value_unitig_pheno_file)
    if not contains_exists or not value_exists:
        # instantiate local vars to be passed to worker processes
        sim = load_pickle(sim_file)
        pim = load_pickle(pim_file)
        q = Manager().Queue()
        # instantiate worker processes to process large unitig file
        process_file(unitig_db, unitig_sample_map_file, sim=sim, pim=pim,
            uim_file=uim_file, q=q, truths=truths_dict)

        # drain queue and write to output files sequentially
        while not q.empty():
            if params.get('truths'):
                unitig_sample_chunk, unitig_pheno_chunk, truths_chunk = q.get()
            else:
                unitig_sample_chunk, unitig_pheno_chunk = q.get()
            if not contains_exists:
                write_list(unitig_sample_chunk, contains_sample_unitig_file)
            if not value_exists:
                write_list(unitig_pheno_chunk, value_unitig_pheno_file)
            if params.get('truths') and not file_exists(truths_outfile):
                unitig_pheno_chunk = [f'{unitig_pheno_chunk[i]}\t{truths_chunk[i]}' for i in range(len(unitig_pheno_chunk))]
                write_list(unitig_pheno_chunk, truths_outfile)

if __name__ == '__main__':
    parse_args()
    main()

