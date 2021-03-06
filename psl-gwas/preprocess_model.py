#!/usr/bin/env python3
###############################################################################
##  preprocess_model.py
##  This file holds helper functions for preprocess.py.
###############################################################################
from random import Random, randint, random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from utility import printd, write_files, complement, file_exists
from collections import Counter

# check if a kmer occurs in fewer than upper samples and
# in more than lower samples. This filters out kmers that are either
# shared by all samples, in which case they cannot cause the pheno,
# or are present in only a few samples, in which case they will not have
# enough statistical power.
def kmer_frequency_fails(val, upper, lower):
    if lower <= int(val) < upper:
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

# consolidates kmers into unitigs. Note that we use python's dict property
# that remembers the order in which items are added. Because of this, kmers
# that occur next to each other are likely to be adjacent, barring any that
# have been filtered out or were broken up by chunkifying. This step is
# solely for reducing the number of variables, so if we miss a few, ie there
# are some kmers that could be consolidated into a unitig but were not,
# it will not affect the accuracy of the computation.
def consolidate(data, k):
    printd('Consolidating chunk...')
    prev_line = data.pop()
    prev_unitig = prev_line[0]
    unitigs = []
    while(data): 
        line = data.pop()
        this_unitig = line[0]

        # kmers are sequential and the same set of samples contain both kmers
        if prev_unitig[0:k - 1] == this_unitig[-(k - 1):] \
                and len(line) == len(prev_line) \
                and set(line[1:]) == set(prev_line[1:]):
            this_unitig = this_unitig[0] + prev_unitig
            line = (this_unitig, *line[1:])
        else:
            unitigs.append(prev_line)
        prev_line = line
        prev_unitig = this_unitig
    printd('Finished consolidating chunk.')
    return unitigs

def format_tuple(tup):
    if type(tup) is str:
        return tup
    else:
        return f'{tup[0]},{tup[1]}'

# removes all unitigs that are not correlated with any pheno.
# thresh is a user-defined value that sets the number of samples which
# contain this unitig and have recorded pheno values.
# prop is the minimum ratio, for any pheno, of samples that do not display
# the pheno to samples that display the pheno, for this unitig to be kept.
def filter_kmers(data, thresh, dfdisp, dfnodisp, kmer_sample_file,
        kmer_pheno_file, lock):
    printd('Filtering kmers...')
    nphenos = dfdisp.shape[1]
    kmer_samples = []
    kmer_phenos = []
    while(data):
        line = data.pop()
        kmer = line[0]
        # collect resistant/vulnerable frequencies for each antibiotic for
        # this unitig
        disp = sum(dfdisp[sample_id[0]] for sample_id in line[1:])
        nodisp = sum(dfnodisp[sample_id[0]] for sample_id in line[1:])

        # 1 test per antibiotic; unitig needs to pass only 1 to avoid
        # getting filtered out
        samples_thresh = 5
        a = np.where((disp + nodisp >= samples_thresh) \
                    & (disp / (disp + nodisp + .01) > thresh))[0]
        if a.size == 0:
            continue
        kmer_pheno_chunk = [kmer]
        for pheno in a:
            kmer_pheno_chunk.append(str(pheno))
        kmer_phenos.append('\t'.join(kmer_pheno_chunk))
        
        kmer_samples.append('\t'.join(map(format_tuple, line)))
        # write every 500K kmers to keep memory consumption under control
        if len(kmer_phenos) >= 500000:
            write_files(lock,
                (kmer_samples, kmer_sample_file),
                (kmer_phenos, kmer_pheno_file))
            kmer_samples = []
            kmer_phenos = []
    write_files(lock,
        (kmer_samples, kmer_sample_file),
        (kmer_phenos, kmer_pheno_file))
    printd('Finished filtering kmers.')
    return
    

# takes a random sample of kmers and creates a distance matrix between
# samples. Optionally takes in a random seed.
def sample_kmers(data, n, seed=randint(1,100000)):
    printd('Sampling kmers...')
    sample_matrix = np.zeros((n, n)) 
    rng = Random(seed)
    num_kmers = int(len(data) * 0.05)
    sampled = rng.sample(data, num_kmers)

    for line in sampled:
        samplelist = line[1:]
        for i, s1 in enumerate(samplelist):
            for s2 in samplelist[i:]:
                sample_matrix[s1[0]][s2[0]] += 1
                sample_matrix[s2[0]][s1[0]] += 1
    printd('Finished sampling kmers.')
    return num_kmers, sample_matrix

# Creates a dictionary of all kmers passed in data, and their complements
# Then, iterates through genomes. If a genome contains a kmer, that genome's
# metadata is added to dict entry for that kmer.
# calls sample_kmers, consolidate, and filter_unitigs
def create_kmer_sample_map(data, raw, q, k, upper, lower, thresh,
        dfdisp, dfnodisp, sim, n, lock, kmer_sample_file, kmer_pheno_file):
    printd('Creating kmer sample map...') 
    # get all kmers in chunk and complement them
    kmers = {}
    for line in data:
        kmer, count = line.split('\t')
        if kmer_frequency_fails(count, upper, lower):
            continue
        kmers[kmer] = Counter()
    
    # map all kmers in chunk to samples containing them
    for count, (raw_id, seq) in enumerate(raw.items()):
        sample_id = sim.get(raw_id, None)
        if sample_id is None:
            continue
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= k: # ensure this contig is long enough to sample
                for i in range(l - k + 1):
                    kmer = contig[i: i + k]
                    kmerlist = kmers.get(kmer, None)
                    if kmerlist is not None:
                        kmerlist[sample_id] += 1
                    else:
                        complist = kmers.get(complement(kmer), None)
                        if complist is not None:
                            complist[sample_id] += 1
    kmers = [(key, *v.items()) for key,v in kmers.items()]
    printd('Finished creating kmer sample map.')
    num_kmers, sample_matrix = sample_kmers(kmers, n)
    printd('Putting data in queue')
    q.put((num_kmers, sample_matrix))
    printd('Finished putting data in queue')
    if kmer_sample_file is None and kmer_pheno_file is None:
        return
    
    # consolidate() will clear kmers list as it builds unitigs list
    # with net 0 memory gain
    kmers = consolidate(kmers, k)
    # filter_unitigs() will clear unitigs list as it builds new unitigs list
    # with net 0 memory gain
    filter_kmers(kmers, thresh, dfdisp, dfnodisp,
        kmer_sample_file, kmer_pheno_file, lock)
    return

## combine similarity matrices from each thread into a single matrix,
## and get total number of kmers sampled
def similar_sample(sample_matrix, num_kmers, similarities_tsv,
        hist_orig_file, hist_sim_scaled_file, hist_dissim_scaled_file,
        similarities_file, dissimilarities_file):
    if not file_exists(similarities_tsv):
        # scale similarities matrix by the mean num sampled kmers each sample
        # shares with itself. Then, normalize to [0,1]. Then remove the diagonal
        # and the lower triangle of the array (since it is symmetric about the
        # major diagonal), and finally round values to 4 decimal places.
        mean_shared_w_self = sample_matrix.diagonal().mean()
        sample_matrix /= mean_shared_w_self
        sample_matrix += 0.001 # ensure all values are nonzero
        sample_matrix *= 1.0/sample_matrix.max()
        np.fill_diagonal(sample_matrix, np.nan)
        sample_matrix = np.triu(sample_matrix)
        np.round(sample_matrix, 4)

        df = pd.DataFrame(sample_matrix)

        # dump to tsv file for ease of restoring, and because tsv file of similarities
        # is a common input to other mGWAS programs
        df.to_csv(similarities_tsv, sep='\t')
    
    else:
        df = pd.read_csv(similarities_tsv, sep='\t', index_col=0)
   
    # create similarity histogram and save it
    plt.hist(df.values, facecolor='green')
    plt.savefig(hist_orig_file, dpi=150)
    plt.clf()
    df = df.stack()
    df = df.reset_index()
    df = df[df[0] > 0] # remove the lower half of the triangle
    # set threshold; 0.75 means drop lowest 75%, keep highest 25%
    highthresh = 0.9
    lowthresh = 0.1
    # find numeric cutoff; the lowest 75% of the data are below this value
    highcutoff = df[0].quantile(highthresh)
    lowcutoff = df[0].quantile(lowthresh)
    # cut off all everything in the middle; only keep the very similar and very dissimilar
    simdf = df[df[0] >= highcutoff].copy(deep=True)
    dissimdf = df[df[0] <= lowcutoff].copy(deep=True)
    dissimdf[0] = 1 - dissimdf[0]
    dfs = (simdf, dissimdf)
    # determine new min, max, range
    files = ((hist_sim_scaled_file, similarities_file), (hist_dissim_scaled_file, dissimilarities_file))
    for i, (pngfile, outfile) in enumerate(files):
        df = dfs[i]
        min_ = df[0].min()
        max_ = df[0].max()
        range_ = max_ - min_
        # shift df left by the min so the new min is 0
        df[0] -= min_
        # rescale data to [0,0.5] or [0,1]
        if i == 0: # high
            scale_factor = 2
            intercept = 0.5
        else: # low
            scale_factor = 2
            intercept = 0.5
        df[0] /= range_ * scale_factor
        # rescale data to [0.5, 1] or [0, 1]
        df[0] += intercept
        # create similarity histogram and save it
        try:
            plt.hist(df[0], bins=50, facecolor='green')
            plt.savefig(pngfile, dpi=150)
            plt.clf()
        except ValueError as e:
            printd('Unable to generate histogram of scaled data')

        # write to csv 
        df.to_csv(outfile, sep='\t', index=False, header=False)

def convert_to_binary(x):
    if x > 1.0:
        return 1
    if x < 0.0:
        return 0
    return x

# separate phenos df into 2 dfs, one holding phenos present and 1 holding
# phenos absent.
def create_disp_nodisp_dfs(phenos, sim):
    df = pd.read_csv(phenos, delimiter='\t')
    idcol = df.columns[0]
    df[idcol] = df[idcol].map(sim)
    df = df[df[idcol].notna()]
    df.set_index(idcol, inplace=True)
    df = df.applymap(convert_to_binary)
    dfdisp = df.fillna(0)
    dfnodisp = 1 - df.fillna(1)
    return dfdisp.to_numpy(), dfnodisp.to_numpy()

