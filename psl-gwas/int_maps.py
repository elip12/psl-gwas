#!/usr/bin/env/python3
###############################################################################
##  int_maps.py
##  This file holds helper functions for creating dict mappings between
##  samples/phenos/kmers and integers. Useful for speeding up PSL.
###############################################################################
import pickle
from utility import printd
import pandas as pd

def create_sample_int_map(samples, phenos, sim_file):
    printd('Creating sample int map...')
    sim = {}
    with open(samples, 'r') as f:
        lines = f.readlines()
    phenosdf = pd.read_csv(phenos, sep='\t', index_col=0)
    phenosdf.dropna(how='all', inplace=True)
    droppedsamples = []
    i = 0
    for line in lines[1:]: # ignore header
        name = line.split('\t')[0]
        if name in phenosdf.index:
            sim[name] = i
            sim[i] = name
            i += 1
        else:
            droppedsamples.append(name)
    with open(sim_file, 'wb') as f:
        pickle.dump(sim, f)
    if len(droppedsamples) > 0:
        printd(('Ignoring samples not present in both'
            f' sample and pheno files: {droppedsamples}'))
    printd('Successfully created sample int map.')

def create_pheno_int_map(phenos, pim_file):
    printd('Creating pheno int map...')
    with open(phenos, 'r') as f:
        line = f.readline() # only need first line which contains headers
    pim = {}
    for i, p in enumerate(line.split()[1:]): # skip id column
        pim[i] = p
        pim[p] = i
    with open(pim_file, 'wb') as f:
        pickle.dump(pim, f)
    printd('Successfully created pheno int map.')

def create_kmer_int_map(kmers, kim_file):
    printd('Creating kmer int map...')
    kim = {}
    with open(kmers, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        u = line.split('\t')[0]
        kim[u] = i
        kim[i] = u
    with open(kim_file, 'wb') as f:
        pickle.dump(kim, f)
    printd('Successfully created kmer int map.')

