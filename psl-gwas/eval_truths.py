import pandas as pd
import numpy as np
from utility import load_pickle, dump_pickle
import random

base = 'ecoli/data/preprocessed/'
#truths_file = base + 'truth_unitig_pheno.txt'
unitig_int_map = load_pickle(base + 'unitig_int_map.pkl')
unitig_sample_map_file = base + 'unitig_sample_map.txt'
phenos_df = pd.read_csv('ecoli/data/raw/phenos.tsv', sep='\t', index_col=0)
sim = load_pickle(base + 'sample_int_map.pkl')
pim = load_pickle(base + 'pheno_int_map.pkl')

def get_truth_kmer_seqs_that_occur_in_data():
    with open(truths_file, 'r') as f:
        truthslines = f.readlines()
    #with open(unique_kmers_file, 'r') as f:
    #    kmerslines = f.readlines()
    #kmerslines = set(k.split('\t')[0] for k in kmerslines)
    truthslines = [l for l in truthslines if l != '\n']
    truthslines = [l.split('\t') for l in truthslines]
    truthslines = [(unitig_int_map[int(k)], p) for k, p,_ in truthslines]
    return truthslines

def read_unitig_sample_map_into_dict():
    with open(unitig_sample_map_file, 'r') as f:
        lines = f.readlines()
    lines = {l.split('\t')[0]:l.split('\t')[1:] for l in lines}
    return lines

def cross_truth_seqs_with_usm(usm, truthslines):
    return [(k,p, usm[k]) for k, p in truthslines]

def convert_sample_ids_to_pheno_values(data):
    means = []
    counts = []
    ps = []
    smps = []
    lengths = []
    for k, p, samples in data:
        if random.random() > 0.01:
           continue
        lengths.append(len(k))
        samplelist = [phenos_df[pim[int(p)]][sim[int(s.split(',')[0])]] for s in samples]
        samplelist = [s for s in samplelist if str(s) != 'nan']
        if len(samplelist) < 1:
            continue
        mean = np.mean(samplelist)
        counts.append(round(mean, 3))
        means.append(round(len(samplelist) / (len(sim)/2), 2))
        ps.append(pim[int(p)])
        #smps.append([sim[int(s.split(',')[0])] for s in samples])
    means = np.array(means)
    
    print('min: ', min(means))
    print('mean:', np.mean(means))
    print('max: ', max(means))
    for i in range(2, 99, 2):
        print(f'{i}%:  ', np.percentile(means, i))
    
    print('mean length:', np.mean(lengths))
    #meancounts = zip(means, counts)
    #meancounts = sorted(meancounts, key=lambda x: x[0])
    #print(meancounts[:30])

def main():
    #truths = get_truth_kmer_seqs_that_occur_in_data()
    #usm = read_unitig_sample_map_into_dict()
    #crossed = cross_truth_seqs_with_usm(usm, truths)
    #dump_pickle(crossed, 'crossed_truth.pkl')
    crossed = load_pickle('crossed_target.pkl')
    convert_sample_ids_to_pheno_values(crossed)

main()
