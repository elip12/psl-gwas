#!/usr/bin/env python3 
from sys import argv

def complement(s):
    m = str.maketrans('ATCG', 'TAGC')
    return s.translate(m)

def create_truths_dict(truths_infile):
    truths = {}
    with open(truths_infile, 'r') as f:
        lines = f.readlines()
    phenos = None
    pim_code = None
    for line in lines:
        if line.startswith('>'):
            phenos = line.rstrip().split('_')
            name = phenos[0]
            phenos = phenos[1:]
            for pheno in phenos:
                if pheno not in truths:
                    truths[pheno] = []
        elif phenos is not None:
            for pheno in phenos:
                truths[pheno].append((name, line.rstrip()))
    return truths

def chunkify(unitig):
    size = 60
    rem = len(unitig) % size
    for i in range(0, len(unitig) - rem, size):
        # if next chunk would be smaller than size, add it to this chunk
        if len(unitig) - (i + size) < size:
            yield unitig[i:]
        else:
            yield unitig[i: i + size]

def unitig_in_truths(unitig, truths, pheno):
    seqs = truths.get(pheno, None)
    if seqs is None:
        return False
    for name, seq in seqs:
        if len(unitig) > 60:
            for unitigchunk in chunkify(unitig):
                comp = complement(unitigchunk)
                if unitigchunk in seq or seq in unitigchunk or comp in seq or seq in comp:
                    return True
        else:
            comp = complement(unitig)
            if unitig in seq or seq in unitig or comp in seq or seq in comp:
                return True
    return False


def get_kmers(pyseer_output, n):
    kmers = []
    with open(pyseer_output, 'r') as f:
        lines = f.readlines()
    lines = lines[:n]
    for line in lines:
        for word in line.split():#','):
            if set(word.rstrip()) == set('ATCG'):
                kmers.append(word.rstrip())
    return kmers
        
def kmer_in_gene(kmer, gene):
    if len(kmer) > 60:
        for kmerchunk in chunkify(kmer):
            comp = complement(kmerchunk)
            if kmerchunk in gene or gene in kmerchunk or comp in gene or gene in comp:
                return True
    else:
        comp = complement(kmer)
        if kmer in gene or gene in kmer or comp in gene or gene in comp:
            return True
    return False

def compute_recall(kmers, truths, pheno):
    genes_covered = 0
    for name, gene in truths[pheno]:
        for kmer in kmers:
            if kmer_in_gene(kmer, gene):
                genes_covered += 1
                print(name)
                break
    rec = genes_covered / len(truths[pheno])
    return round(rec, 5)

def compute_precision(kmers, truths, pheno):
    tps = 0
    for kmer in kmers:
        if unitig_in_truths(kmer, truths, pheno):
            tps += 1
    pre = tps / len(kmers)
    return round(pre, 5)

def compute_rank_of_1st_hit(kmers, truths, pheno):
    for i, kmer in enumerate(kmers): # guarantee to be descending order
        if unitig_in_truths(kmer, truths, pheno):
            return round(1 / (i + 1), 5)

def main():
    truths_infile = argv[1]
    pyseer_output = argv[2]
    pheno = pyseer_output.split('.')[1].capitalize()
    truths = create_truths_dict(truths_infile)
    kmers = get_kmers(pyseer_output, int(argv[3]))
    error = compute_precision(kmers, truths, pheno)
    print(f'Precision ({pheno}) (kmers):\t{error}')
    genes_missed = compute_recall(kmers, truths, pheno)
    print(f'Hits ({pheno}) (genes):\t{genes_missed}')
    rank = compute_rank_of_1st_hit(kmers, truths, pheno)
    print(f'MRR ({pheno}):\t{rank}')

if __name__ == '__main__':
    main()
