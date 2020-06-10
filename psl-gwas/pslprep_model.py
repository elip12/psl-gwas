#!/usr/bin/env python3
###############################################################################
##  pslprep_model.py
##  This file holds helper functions for pslprep.py
###############################################################################
from utility import write_files, write_list, load_pickle, printd, complement
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
    score = None
    for line in lines:
        if line.startswith('>'):
            linelist = line.rstrip().split('_')
            if 'baseline' in linelist[0]:
                baseline = linelist[0].split('-')
                if len(baseline) == 2:
                    try:
                        score = float(baseline[1])
                        if score < 0.0 or score > 1.0:
                            score = None
                    except ValueError as e:
                        score = None
            phenos = linelist[1:]
            for pheno in phenos:
                pim_code = pim.get(pheno, None)
                if pim_code is None:
                    print('Error: Pheno in truths does not exist in \
                    phenos file:', pheno)
                    print('Exiting')
                    exit(1)
                if pim_code not in truths:
                    truths[pim_code] = []
        elif phenos is not None:
            for pheno in phenos:
                pim_code = pim[pheno]
                truths[pim_code].append((line.rstrip(), score))
    return truths
        
def kmer_in_truths(kmer, truths, pheno):
    seqs = truths.get(int(pheno), None)
    if seqs is None:
        return False
    comp = complement(kmer)
    for seq, score in seqs:
        if kmer in seq or seq in kmer or comp in seq or seq in comp:
            if score is None:
                return True
            else:
                return score
    return False
    
def pim_truths(pim, truths):
    for pheno in truths:
        if pheno in pim:
            yield pim[pheno], pheno

def kmer_pheno_db(data, kim_file, value_kmer_pheno_file,
        truth_kmer_pheno_file, baseline_kmer_pheno_file, lock, truths=None,
        baseline=None):
    kim = load_pickle(kim_file)
    kmer_pheno_chunk = []
    if truths:
        truths_chunk = []
    else:
        truths_chunk = None
    if baseline:
        baseline_chunk = []
    else:
        baseline_chunk = None
    for line in data:
        linelist = line.split('\t')
        kmer = linelist[0]
        for pheno in linelist[1:]:
            kmer_pheno_chunk.append(f'{kim[kmer]}\t{pheno}')
            if truths and kmer_in_truths(kmer, truths, pheno):
                truths_chunk.append(f'{kim[kmer]}\t{pheno}')
            if baseline:
                score = kmer_in_truths(kmer, baseline, pheno)
                if score is True:
                    baseline_chunk.append(f'{kim[kmer]}\t{pheno}')
                elif score is not False and score > 0.0 and score < 1.0:
                    baseline_chunk.append(f'{kim[kmer]}\t{pheno}\t{score}')
        if len(kmer_pheno_chunk) >= 500000:
            write_files(lock,
                (kmer_pheno_chunk, value_kmer_pheno_file),
                (truths_chunk, truth_kmer_pheno_file),
                (baseline_chunk, baseline_kmer_pheno_file))
            kmer_pheno_chunk = []
            if truths_chunk is not None:
                truths_chunk = []
    write_files(lock,
        (kmer_pheno_chunk, value_kmer_pheno_file),
        (truths_chunk, truth_kmer_pheno_file),
        (baseline_chunk, baseline_kmer_pheno_file))

# convert kmer db to psl input
def kmer_sample_db(data, kim_file, contains_sample_kmer_file, lock, truths=None):
    kim = load_pickle(kim_file)
    kmer_sample_chunk = []
    for line in data:
        linelist = line.split('\t')
        kmer = linelist[0]
        
        kmer_sample_lines = []
        for sample_ in linelist[1:]:
            sample = sample_.split(',')
            kmer_sample_lines.append(f'{sample[0]}\t{kim[kmer]}') # add some perturbation of sample[1] to the end to get different truth values for different num CNVs
        kmer_sample_chunk.append('\n'.join(kmer_sample_lines))
        
        # write every 500k to limit memory usage
        if len(kmer_sample_chunk) >= 500000:
            lock.acquire()
            write_list(kmer_sample_chunk, contains_sample_kmer_file)
            lock.release()
            kmer_sample_chunk = []
    write_files(lock, (kmer_sample_chunk, contains_sample_kmer_file))

# convert phenos tsv to psl input
def sample_pheno(phenos, sim, pim, outfile):
    df = pd.read_csv(phenos, delimiter='\t')
    
    idcol = df.columns[0]
    df['sample_id'] = df[idcol].apply(lambda x: sim.get(x, -1))    # int(sim.get(x)) if sim.get(x) is not None else np.nan)
    df = df[df['sample_id'] != -1]
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
    df = df[df['value'] > 0]
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
