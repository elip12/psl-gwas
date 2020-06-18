#!/usr/bin/env python3
###############################################################################
##  postprocess_model.py
##  This file provides helper functions for the postprocessing stage. 
###############################################################################
from os.path import basename, join
import pandas as pd
from utility import printd, write_list

def process(data, lock, pim, kim_file, thresh, fsa_file, scored_kmers_file):
    kim = load_pickle(kim_file)
    chunk = []
    for line in data:
        linelist = line.split()
        if float(linelist[2]) < thresh:
            continue
        outline = (kim[int(linelist[0])], pim[int(linelist[1])], linelist[2])
        chunk.append(outline)
    kmers = [f'>{i}\n{line[0]}' for i, line in enumerate(chunk)]
    values = ['\t'.join(tup) for tup in chunk]
    write_files(lock,
        (values, scored_kmers_file),
        (kmers, fsa_file))

def consolidate_model(data, k, p=False):
    for i, (unitig1, ranks1) in enumerate(data):
        for j in range(i + 1, len(data)):
            unitig2, ranks2 = data[j]
            if unitig1[:k] == unitig2[-k:]:
                newunitig = unitig2[:-k] + unitig1
                if p:
                    print(j, '\t', unitig2)
                    print(i, '\t', unitig1.rjust(len(newunitig), ' '))
                    print('-\t', '-'.rjust(len(newunitig), '-'))
                    print(j, '\t', newunitig, '\n')
                data[i] = None
                data[j] = (newunitig, ranks1 + ranks2)
                break
            elif unitig1[-k:] == unitig2[:k]:
                newunitig = unitig1[:-k] + unitig2
                if p:
                    print(i, '\t', unitig1)
                    print(j, '\t', unitig2.rjust(len(newunitig), ' '))
                    print('-\t', '-'.rjust(len(newunitig), '-'))
                    print(j, '\t', newunitig, '\n')
                data[j] = (newunitig, ranks1 + ranks2)
                data[i] = None
                break
    return [d for d in data if d is not None]


def consolidate(name, unitigs, outdir):
    unitigs = list(unitigs.iloc[:, 0])
    unitigs = [(u, [1/(i+1)]) for i, u in enumerate(unitigs)]
    printd('Original num kmers:', len(unitigs))
    for k in range(30, 20, -1):
        unitigs = consolidate_model(unitigs, k, p=False)
        printd(f'K: {k}, num kmers: {len(unitigs)}')
    printd('Final num kmers:', len(unitigs))
    unitigs = [(unitig, 1 / (sum(ranks) / len(ranks))) for unitig, ranks in unitigs]
    unitigs = sorted(unitigs, key=lambda x: x[1])
    unitigs = [u[0] for u in unitigs]
    write_list(unitigs, join(outdir, name))

def separate_phenos(infile, outdir, n, no_consolidate=False)
df = pd.read_csv(infile, sep='\t', header=None)
for col in df[1].unique():
    dfcol = df[df[1] == col]
    dfcol = dfcol.nlargest(n, columns=2)
    name = f'psl.{col.lower()}'
    if no_consolidate:
        dfcol.to_csv(join(outdir, name), sep='\t')
    else:
        consolidate(name, dfcol, outdir)
