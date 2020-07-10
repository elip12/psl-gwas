#!/usr/bin/env python3 
import pandas as pd
import numpy as np
from sys import argv
from os import listdir

def normalize(pvals):
    coef = 0.1
    pvals = np.array(pvals)
    pvals = 1 - pvals
    range_ = 1 - pvals.min()
    if range_ != 0:
        pvals -= pvals.min()
        pvals /= range_ * (1 / coef)
        pvals += 1 - coef
    pvals = np.around(pvals, 4)
    return list(pvals)

for fname in listdir(argv[1]):
    if fname.startswith('pyseer'):
        df = pd.read_csv(f'{argv[1]}/{fname}', sep='\t')
        col = fname.split('.')[1]
        df = df.nsmallest(int(argv[3]), columns='lrt-pvalue')
        variants = list(df['variant'])
        pvals = list(df['lrt-pvalue'])
        pvals = normalize(pvals)
        variants = list(zip(variants, pvals))
        variants = (f'>baseline-{score}_{col.capitalize()}\n{unitig}\n' for unitig, score in variants)
        with open(f'{argv[2]}/baseline.fa', 'a+') as f:
            f.writelines(variants)
