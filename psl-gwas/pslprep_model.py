#!/usr/bin/env python3
###############################################################################
##  pslprep_model.py
##  This file holds helper functions for pslprep.py
###############################################################################
from utility import write_2_files, write_list, load_pickle, printd
import pandas as pd
import numpy as np
import pickle
from sys import argv

def create_truths_dict(truths_infile, pim):
    truths = {}
    with open(truths_infile, 'r') as f:
        lines = f.readlines()
    phenos = None
    pim_code = None
    for line in lines:
        if line.startswith['>']:
            phenos = line.rstrip().split('_')[1:]
            for pheno in phenos:
                pim_code = pim.get(pheno, None)
                if pim_code is None:
                    print('Error: Pheno in truths does not exist in phenos file')
                    print('Exiting')
                    exit(1)
                if pim_code not in truths:
                    truths[pim_code] = []
        elif phenos is not None and pim_code is not None:
            for pheno in phenos:
                truths[pim_code].append(line)
    return truths
        
def unitig_in_truths(unitig, truths, pheno):
    seqs = truths.get(pheno, None)
    if seqs is None:
        return False
    for seq in seqs:
        if unitig in seq or seq in unitig:
            return 1
    return 0
    
def pim_truths(pim, truths):
    for pheno in truths:
        if pheno in pim:
            yield pim[pheno], pheno

def unitig_pheno_db(data, uim_file, value_unitig_pheno_file,
        truth_unitig_pheno_file, lock, truths=None):
    uim = load_pickle(uim_file)
    unitig_pheno_chunk = []
    if truths:
        truths_chunk = []
    else:
        truths_chunk = None
    for line in data:
        linelist = line.split('\t')
        unitig = linelist[0]

        for pheno in linelist[1:]:
            unitig_pheno_chunk.append(f'{uim[unitig]}\t{pheno}')
            if truths:
                if unitig_in_truths(unitig, truths, pheno):
                    truths_chunk.append(f'{uim[unitig]}\t{pheno}\t1')
        if len(unitig_pheno_chunk) >= 500000:
            write_2_files(unitig_pheno_chunk, value_unitig_pheno_file,
                truths_chunk, truth_unitig_pheno_file, lock)
            unitig_pheno_chunk = []
            truths_chunk = []
    write_2_files(unitig_pheno_chunk, value_unitig_pheno_file,
        truths_chunk, truth_unitig_pheno_file, lock)

# convert unitig db to psl input
def unitig_sample_db(data, uim_file, contains_sample_unitig_file, lock, truths=None):
    uim = load_pickle(uim_file)
    unitig_sample_chunk = []
    for line in data:
        linelist = line.split('\t')
        unitig = linelist[0]
        
        unitig_sample_lines = []
        for sample_ in linelist[1:]:
            sample = sample_.split(',')
            unitig_sample_lines.append(f'{sample[0]}\t{uim[unitig]}\tsample[1]')
        unitig_sample_chunk.append('\n'.join(unitig_sample_lines))
        
        # write every 500k to limit memory usage
        if len(unitig_sample_chunk) >= 500000:
            lock.acquire()
            write_list(unitig_sample_chunk, contains_sample_unitig_file)
            lock.release()
            unitig_sample_chunk = []
    lock.acquire()
    write_list(unitig_sample_chunk, contains_sample_unitig_file)
    lock.release()

# convert phenos tsv to psl input
def sample_pheno(phenos, sim, pim, outfile):
    df = pd.read_csv(phenos, delimiter='\t')
    
    idcol = df.columns[0]
    df['sample_id'] = df[idcol].apply(lambda x: sim[x])    # int(sim.get(x)) if sim.get(x) is not None else np.nan)
    #df['sample_id'] = df['sample_id'].dropna()
    df.drop(idcol, axis=1, inplace=True)
    df.rename(columns=pim, inplace=True) 
    
    new = pd.DataFrame(columns=['sample_id', 'pheno', 'value'])
    for col in df:
        if col != 'sample_id':
            slice_ = df[['sample_id', col]].copy()
            slice_['pheno'] = col
            slice_.rename(columns={col: 'value'}, inplace=True)
            new = pd.concat([new, slice_[['sample_id', 'value', 'pheno']]],
                ignore_index=True, sort=False) 
    df = new.dropna()
    df.to_csv(outfile, sep='\t', index=False, header=False)

# rescale similar pheno data to have better distribution
# cut off everything below 0.75, and stretch the remaining 25%
# to be 0.5 <= x <= 1
def remap_similar_pheno(series):
    return series.apply(lambda x: 0 if x < 0.75 else 1 - ((1 - x) * 2))

# find correlations between phenos from input data and create pheno similarity
# file 
def similar_pheno(phenos, pim, outfile):
    df = pd.read_csv(phenos, delimiter='\t')
    
    idcol = df.columns[0]
    df.drop(idcol, axis=1, inplace=True)  
    df.rename(columns=pim, inplace=True) 
    df = df.corr()
    df = df.stack()
    df = df.reset_index()
    # turn everything under .75 into 0, and remap everything remaining to [0.5,1]
    df[0] = remap_similar_pheno(df[0])
    df.to_csv(outfile, sep='\t', index=False, header=False)


if __name__ == '__main__':
    parse_args()
    main_wrapper()
