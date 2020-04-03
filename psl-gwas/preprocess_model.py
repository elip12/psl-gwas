#!/usr/bin/env python3
###############################################################################
##  preprocess_model.py
##  This file holds helper functions for preprocess.py.
###############################################################################
from random import Random, randint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# global complement map, only needs to be created once
COMPLEMENT_MAP = str.maketrans('ATCG', 'TAGC')

# check if a kmer occurs in fewer than upper samples and
# in more than lower samples. This filters out kmers that are either
# shared by all samples, in which case they cannot cause the pheno,
# or are present in only a few samples, in which case they will not have
# enough statistical power.
def kmer_frequency_fails(val, upper, lower):
    if lower <= int(val) <= upper:
        return False
    return True

# extracts the number of samples in the GWAS from the samples file
def num_samples(fname):
    with open(fname, 'r') as f:                                                
        return len(f.readlines()) - 1

# reads in all contigs from all input files into a dictionary
def parse_input(infile):
    seqs = {}
    with open(infile, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        sample, fname = line.rstrip().split('\t')
        if not fname.endswith('.fa'):
            continue
        with open(fname, 'r') as f:
            contiglines = f.readlines()
        seqs[sample] = [contigline.rstrip() for contigline in contiglines \
            if contigline[0] in ['A','T','C','G']]
    return seqs 

# complements a kmer
def complement(kmer):
   return kmer.translate(COMPLEMENT_MAP)

# consolidates kmers into unitigs. Note that we use python's dict property
# that remembers the order in which items are added. Because of this, kmers
# that occur next to each other are likely to be adjacent, barring any that
# have been filtered out or were broken up by chunkifying. This step is
# solely for reducing the number of variables, so if we miss a few, ie there
# are some kmers that could be consolidated into a unitig but were not,
# it will not affect the accuracy of the computation.
def consolidate(data, k):
    prev_line = data[0]
    prev_unitig = prev_line[0]
    prev_samples = prev_line[1]
    unitig_db_chunk = []
    for line in data[1:]:
        this_unitig = line[0]
        this_samples = line[1]

        # kmers are sequential and the same set of samples contain both kmers
        if prev_unitig[-(k - 1):] == this_unitig[0:k - 1] \
                and len(this_samples) == len(prev_samples) \
                and set(this_samples) == set(prev_samples):
            this_unitig = f'{prev_unitig}{this_unitig[-1]}'
            prev_line = (this_unitig, prev_line[1])
        else:
            unitig_db_chunk.append(prev_line)
            prev_line = line
        
        prev_unitig = this_unitig
        prev_samples = this_samples
    return unitig_db_chunk

# removes all unitigs that are not correlated with any pheno.
# thresh is a user-defined value that sets the number of samples which
# contain this unitig and have recorded pheno values.
# prop is the minimum ratio, for any pheno, of samples that do not display
# the pheno to samples that display the pheno, for this unitig to be kept.
def filter_unitigs(data, thresh, dfdisp, dfnodisp, prop=0.05):
    unitig_db_chunk = []
    nphenos = dfdisp.shape[1]
    for linelist in data:
        unitig_db_line = [linelist[0]]
        disp = np.zeros(nphenos)
        nodisp = np.zeros(nphenos)

        for sample_id, _ in linelist[1]:
            if sample_id not in dfdisp.index:
                continue
            unitig_db_line.append(sample_id)
            # collect resistant/vulnerable frequencies for each antibiotic for
            # this unitig
            disp += dfdisp.loc[sample_id].to_numpy()
            nodisp += dfnodisp.loc[sample_id].to_numpy()

        # 1 test per antibiotic; unitig needs to pass only 1 to avoid
        # getting filtered out
        a = np.where((disp + nodisp >= thresh) \
                    & (disp > 0) \
                    & (nodisp / (disp + 0.0001) < prop))[0]
        if a.size > 0:
            unitig_db_chunk.append('\t'.join(unitig_db_line))
    return unitig_db_chunk

# takes a random sample of kmers and creates a distance matrix between
# samples. Optionally takes in a random seed.
def sample_kmers(data, sim, n, seed=None):
    sample_matrix = np.zeros((n, n)) 
    if seed is not None:
        rng = Random(seed)
    else:
        rng = Random(randint(1,100000))

    for line in data:
        if rng.random() > 0.01:
            continue
        for i, s1_ in enumerate(line[:-1]):
            s1 = s1_[0]
            for s2_ in line[i + 1:]:
                s2 = s2_[0]
                sample_matrix[int(sim[s1])][int(sim[s2])] += 1 
                sample_matrix[int(sim[s2])][int(sim[s1])] += 1
    return num_kmers, sample_matrix

# Creates a dictionary of all kmers passed in data, and their complements
# Then, iterates through genomes. If a genome contains a kmer, that genome's
# metadata is added to dict entry for that kmer.
# calls sample_kmers, consolidate, and filter_unitigs
def create_unitig_sample_map(data, raw, k, q, upper, lower, thresh,
        dfdisp, dfnodisp, sim, n):
   
    # get all kmers in chunk and complement them
    kmers = {}
    for line in data:
        kmer, count = line.split('\t')
        if kmer_frequency_fails(count, upper, lower):
            continue
        comp = complement(kmer)
        kmers[kmer] = []
        kmers[comp] = []
    
    # map all kmers in chunk to samples containing them
    for count, (raw_id, seq) in enumerate(raw.items()):
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= k: # ensure this contig is long enough to sample
                for i in range(l - k + 1):
                    kmer = contig[i: i + k]
                    if kmer in kmers:
                        kmers[kmer].append((raw_id, c_id))
    kmers = {k:v for k,v in kmers.items() if len(v) > 0} 
    num_kmers, sample_matrix = sample_kmers(kmers.values(), sim, n)
    unitigs = consolidate(kmers.items(), k)
    unitigs = filter_unitigs(unitigs, thresh, dfdisp, dfnodisp)
    q.put((unitigs, num_kmers, sample_matrix))

    ## combine similarity matrices from each thread into a single matrix,
    ## and get total number of kmers sampled
def similar_sample(sample_matrix, num_kmers, similarities_tsv,
        hist_orig_file, hist_scaled_file, outfile):
    np.fill_diagonal(sample_matrix, np.nan) 
    sample_matrix = np.triu(sample_matrix)

    # scale similarity counts to values in [0, 1]
    sample_matrix /= num_kmers

    df = pd.DataFrame(sample_matrix)

    # dump to tsv file for ease of restoring, and because tsv file of similarities
    # is a common input to other mGWAS programs
    df.to_csv(similarities_tsv, sep='\t')
    
    # optionally read csv for ease of restoring
    # df = pd.read_csv('data/intermediate/similarities.tsv', sep='\t', index_col=0)
   
    # create similarity histogram and save it
    plt.hist(df.values, facecolor='green')
    plt.savefig(hist_orig_file, dpi=150)

    df = df.stack()
    df = df.reset_index()
    
    # set threshold; 0.75 means drop lowest 75%, keep highest 25%
    thresh = 0.9
    # find numeric cutoff; the lowest 75% of the data are below this value
    cutoff = df[0].quantile(thresh)
    # cut off all values below (less similar than) cutoff
    df = df[df[0] > cutoff]
    # determine new min, max, range
    min_ = df[0].min()
    max_ = df[0].max()
    range_ = max_ - min_
    # shift df left by the min so the new min is 0
    df[0] -= min_
    # rescale data to [0,0.5]
    df[0] /= range_ * 2
    # shift right by 0.5 so the new range is [0.5, 1]
    df[0] += 0.5

    # create similarity histogram and save it
    plt.hist(df[0], facecolor='green')
    plt.savefig(hist_scaled_file, dpi=150)

    # write to csv 
    df.to_csv(outfile, sep='\t', index=False, header=False)

def convert_to_binary(x):
    if x > 1.0:
        return 1.0
    if x < 0.0:
        return 0.0
    return x

# separate phenos df into 2 dfs, one holding phenos present and 1 holding
# phenos absent.
def create_disp_nodisp_dfs(phenos):
    df = pd.read_csv(phenos, delimiter='\t')
    idcol = df.columns[0]
    df.set_index(idcol, inplace=True)
    
    df = df.applymap(convert_to_binary)
    dfdisp = df.fillna(0)
    dfnodisp = 1 - df.fillna(1)
    return dfdisp, dfnodisp

