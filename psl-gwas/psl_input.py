#!/usr/bin/env python3
from utility import process_file, write_list, parse_args, load_pickle, printd, \
get_params, check_outfile
from create_psl_input_helpers import unitig_db, sample_pheno, similar_pheno
from multiprocessing import Manager, Queue
from os.path import join

def main():
    # get params
    params = get_params()

    # get project name
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
    
    # check output files do not exist
    check_outfile(contains_sample_unitig_file)
    check_outfile(value_sample_pheno_file)
    check_outfile(value_unitig_pheno_file)
    check_outfile(similar_pheno_pheno_file)
    
    # create smaller psl input files that can be efficiently done w 1 thread
    sample_pheno(phenos_file, sim, pim, value_sample_pheno_file)
    similar_pheno(phenos_file, pim, similar_pheno_pheno_file)
    
    # instantiate local vars to be passed to worker processes
    sim = load_pickle(sim_file)
    pim = load_pickle(pim_file)
    q = Manager().Queue()
    # instantiate worker processes to process large unitig file
    process_file(unitig_db, unitig_sample_map_file, sim=sim, pim=pim, uim_file=uim_file, q=q)

    # drain queue and write to output files sequentially
    while not q.empty():
        unitig_sample_chunk, unitig_pheno_chunk = q.get()
        write_list(unitig_sample_chunk, contains_sample_unitig_file)
        write_list(unitig_pheno_chunk, value_unitig_pheno_file)

if __name__ == '__main__':
    parse_args()
    wrapper()

