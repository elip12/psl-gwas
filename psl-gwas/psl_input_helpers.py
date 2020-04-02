#!/usr/bin/env python3
from utility import write_list, load_pickle, printd,
import pandas as pd
import pickle
from sys import argv

# kmer db to psl
def unitig_db(data, sim, pim, uim_file, q):
    uim = load_pickle(uim_file)
    unitig_sample_chunk = []
    unitig_pheno_chunk = []
    for line in data:
        linelist = line.split()
        unitig = linelist[0]
        
        unitig_sample_lines = []
        for sample in linelist[1:]:
            unitig_sample_lines.append(f'{sim[sample]}\t{kim[kmer]}\t1.0')
        unitig_sample_chunk.append('\n'.join(unitig_sample_lines))
        
        num_phenos = int(len(pim) / 2)
        unitig_pheno_lines = [f'{uim[unitig]}\t{i}' for i in range(num_phenos)]
        unitig_pheno_chunk.append('\n'.join(unitig_pheno_lines))
    q.put((unitig_sample_chunk, unitig_pheno_chunk))

# resistance sample class
# outfile is value_sample_pheno
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

# similar antibiotic
def remap_similar_pheno(series):
    return series.apply(lambda x: 0 if x < 0.75 else 1 - ((1 - x) * 2))

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
