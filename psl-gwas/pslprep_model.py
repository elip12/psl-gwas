#!/usr/bin/env python3
###############################################################################
##  pslprep_model.py
##  This file holds helper functions for pslprep.py
###############################################################################
from utility import write_list, load_pickle, printd
import pandas as pd
import pickle
from sys import argv

def create_truths_dict(truths_infile):
    truths = {}
    with open(truths_infile, 'r') as f:
        lines = f.readlines()
    for i in range(0, len(lines), 2):
        gene = lines[i][1:].rstrip()
        seq = lines[i + 1].rstrip()
        truths[gene] = seq
    return truths
        

def unitig_in_truths(unitig, seq):
    if unitig in seq:
        return 1
    return 0
    
def pim_truths(pim, truths):
    for pheno in truths:
        if pheno in pim:
            yield pim[pheno], pheno

# convert unitig db to psl input
def unitig_db(data, sim, pim, uim_file, truths=None, q):
    uim = load_pickle(uim_file)
    unitig_sample_chunk = []
    unitig_pheno_chunk = []
    if truths:
        truth_chunk = []
    for line in data:
        linelist = line.split()
        unitig = linelist[0]
        
        unitig_sample_lines = []
        for sample in linelist[1:]:
            unitig_sample_lines.append(f'{sim[sample]}\t{uim[unitig]}\t1.0')
        unitig_sample_chunk.append('\n'.join(unitig_sample_lines))
        
        if truths is None:
            num_phenos = int(len(pim) / 2)
            unitig_pheno_lines = [f'{uim[unitig]}\t{i}' for i in range(num_phenos)]
        else:
            unitig_pheno_lines = [(f'{uim[unitig]}\t{i}',
                unitig_in_truths(unitig, truths[pheno])
                ) for i, pheno in pim_truths(pim, truths)]
            truth_values = [l[1] for l in unitig_pheno_lines]
            truth_chunk.extend(truth_values)
            unitig_pheno_lines = [l[0] for l in unitig_pheno_lines]
            unitig_pheno_chunk.extend(unitig_pheno_lines)
    if truths:
        q.put((unitig_sample_chunk, unitig_pheno_chunk, truth_chunk))
    else:
        q.put((unitig_sample_chunk, unitig_pheno_chunk))

# convert phenos tsv to psl input
def sample_pheno(phenos, sim, pim, outfile):
    df = pd.read_csv(phenos, delimiter='\t')
    
    idcol = df.columns[0]
    df['sample_id'] = df[idcol].apply(lambda x: sim[x])
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
