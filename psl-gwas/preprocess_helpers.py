from random import sample
import numpy as np
import pandas as pd

# global complement map, only needs to be created once
COMPLEMENT_MAP = str.maketrans('ATCG', 'TAGC')

def kmer_frequency(val, upper, lower):
    if lower <= x <= upper:
        return False
    return True

# for small files
def num_samples(fname):
    with open(fname, 'r') as f:                                                
        return len(f.readlines()) - 1

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

# DNA base complement
def complement(kmer):
   return kmer.translate(COMPLEMENT_MAP)

def consolidate(data, k):
    prev_linelist = data[0]
    prev_unitig = prev_linelist[0]
    prev_samples = prev_linelist[1]
    unitig_db_chunk = []
    for linelist in data[1:]:
        this_unitig = linelist[0]
        this_samples = linelist[1]
        
        # if kmers are sequential and the same set of samples contain both kmers
        if prev_unitig[-(k - 1):] == line[0:k - 1] \
                and len(this_samples) == len(prev_samples) \
                and set(this_samples) == set(prev_samples):
            this_unitig = f'{prev_unitig}{line[k - 1]}'
            prev_line = (this_unitig, prev_line[1])
        else:
            unitig_db_chunk.append(prev_line)
            prev_line = linelist
        
        prev_unitig = this_unitig
        prev_samples = this_samples
    return unitig_db_chunk

def filter_unitigs(data, thresh, dfdisp, dfnodisp, prop=0.05):
    kmer_db_chunk = []
    nphenos = dfdisp.shape[1]
    for linelist in data:
        kmer_db_line = [linelist[0]]
        disp = np.zeros(nphenos)
        nodisp = np.zeros(nphenos)

        for sample_id, _ in linelist[1]:
            if sample_id not in dfdisp.index:
                continue
            kmer_db_line.append(sample_id)
            # collect resistant/vulnerable frequencies for each antibiotic for
            # this kmer
            disp += dfdisp.loc[sample_id].to_numpy()
            nodisp += dfnodisp.loc[sample_id].to_numpy()

        # 1 test per antibiotic; kmer needs to pass only 1 to avoid
        # getting filtered out
        a = np.where((disp + nodisp >= thresh) \
                    & (disp > 0) \
                    & (nodisp / disp < prop))[0]

        if a.size > 0:
            kmer_db_chunk.append('\t'.join(kmer_db_line))
    return kmer_db_chunk

# data is just values in dict
def sample_kmers(data, sim, n)
    sample_matrix = np.zeros((n, n), dtype=np.uint32)
    num_kmers = int(len(data) * 0.01)
    sampled_kmers = sample(data, num_kmers)

    for line in sampled_kmers:
        for i, s1_ in enumerate(line[:-1]):
            s1 = s1_[0]
            for s2_ in linelist[i + 1:]:
                s2 = s2[0]
                sample_matrix[int(sim[s1])][int(sim[s2])] += 1 
                sample_matrix[int(sim[s2])][int(sim[s1])] += 1
    return num_kmers, sample_matrix

# Creates a dictionary of all kmers passed in data, and their complements
# Then, iterates through genomes. If a genome contains a kmer, that genome's
# metadata is added to dict entry for that kmer.
# After all samples have been iterated over, writes kmers to file.
def create_unitig_sample_map(data, raw, k, q, upper, lower, thresh,
        dfdisp, dfnodisp, sim, n):
    # get all kmers in chunk and complement them
    kmers = {}
    for line in data:
        kmer, count = line.split('\t')
        if kmer_frequency_fails(count, upper, lower):
            continue
        comp = complement_kmer(kmer)
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
    num_kmers, sample_matrix = sample_kmers(kmers.values(), sim, n)
    unitigs = consolidate(kmers.items(), k)
    unitigs = filter_unitigs(unitigs, thresh, dfdisp, dfnodisp)
    q.put(unitigs, num_kmers, sample_matrix)

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
    df.to_csv(similarities_file, sep='\t')
    
    # optionally read csv for ease of restoring
    # df = pd.read_csv('data/intermediate/similarities.tsv', sep='\t', index_col=0)
   
    # create similarity histogram and save it
    plt.hist(df.values, bins=10, facecolor='green')
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
    plt.hist(df[0], bins=10, facecolor='green')
    plt.savefig(hist_scaled_file, dpi=150)

    # write to csv 
    df.to_csv(outfile, sep='\t', index=False, header=False)

