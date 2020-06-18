#!/usr/bin/env python3
###############################################################################
##  postprocess.py
##  This script separates the scored kmers into one file per phenotype,
##  then optionally consolidates the best kmers.
###############################################################################
from utility import process_file, write_list, parse_args, load_pickle, \
get_params, file_exists, write_files
from multiprocessing import Manager
from os.path import join
from postprocess_model import process, separate_phenos


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
    outdir = join(project, 'data', 'postprocessed')

    # create output files if they do not exist
    if file_exists(fsa_file):
        fsa_file = None
    if file_exists(scored_kmers_file):
        scored_kmers_file = None
    if fsa_file or scored_kmers_file:
        lock = Manager().Lock()
        pim = load_pickle(pim_file)
        
        process_file(process, INPUT_FILE, lock=lock, pim=pim, uim_file=uim_file,
                fsa_file=fsa_file, scored_unitigs_file=scored_unitigs_file)
    separate_phenos(scored_kmers_file, outdir,
        params['separate-phenos'], params['no-consolidate'])

if __name__ == '__main__':
    parse_args()
    main()

