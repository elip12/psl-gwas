import pickle
from math import sqrt

# read in from kmer sample map reduced
# for each kmer:
#   compute the mean squared distance between all pairs of samples that have that kmer
#   can do l1, l2, linf. maybe start with linf: the max distance between any two pair of samples
# rank kmers by distance, and include a confidence score that is corresponds to normalized distance

def mean_relatedness(m, samples):
    relatednesses = []
    for s1 in samples:
        for s2 in samples:
            if s1 is not s2:
                relatednesses.append(m[s1][s2])
    return sum(relatednesses) / len(relatednesses)

def min_relatedness(m, samples):
    relatednesses = []
    for s1 in samples:
        for s2 in samples:
            if s1 is not s2:
                relatednesses.append(m[s1][s2]) # ensure we are not double counting here.
                #but if we are its probably fine?
                #i think the matrix is only filled in the upper triangle anyway
    return min(relatednesses)

inpath = 'data/intermediate'
outpath = 'gwas/data'
kmers = []
# fill distance matrix
with open(f'{inpath}/relatedness_matrix.pickle', 'rb') as f:
    m = pickle.load(f)
with open(f'{inpath}/sample_int_map.pickle', 'rb') as f:
    sim = pickle.load(f)
with open(f'{inpath}/kmer_sample_map_reduced.txt', 'r') as f:
    lines = f.readlines()
for line_ in lines:
    line = line_.split()
    kmer = line[0]
    samples = [sim[sample] for sample in line[1:]]
    
    # use the int map also
    kmers.append((kmer, min_relatedness(m, samples))
    #kmers.append((kmer, mean_relatedness(m, samples))

kmers = sorted(kmers, key=lambda k: k[1])
with open(f'{outpath}/scored_kmers.txt', 'w') as f:
    f.writelines(f'{kmer}\t{score}\n' for kmer, score in kmers)

