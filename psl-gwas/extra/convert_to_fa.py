import pandas as pd
import numpy as np
from sys import argv
from os import listdir

def normalize(pvals):
    coef = 0.5
    # take the log so that really small p values are negative with large abs values
    pvals = np.log(pvals)
    # normalize to [0, 1]
    range_ = pvals.max() - pvals.min()
    pvals -= pvals.min()
    pvals /= range_
    # invert values so 1 is good and 0 is bad
    pvals = 1 - pvals
    # scale to [0, coef]
    pvals *= coef
    # shift to [1 - coef, 1]
    pvals += 1 - coef
    # round to 4 decimal places
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
