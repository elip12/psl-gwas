import pickle
from os.path import getsize

def create_kmer_int_map():
    kim = {}
    fname = '../data/kmer_sample_map_reduced.txt'
    size = getsize(fname)
    with open(fname, 'r') as f:
        lines = f.read(size).splitlines()
        for i, line in enumerate(lines):
            k = line.split()[0]
            kim[k] = i
            kim[i] = k
    with open('../data/kmer_int_map.pickle', 'wb+') as f:
        pickle.dump(kim, f)
create_kmer_int_map()
