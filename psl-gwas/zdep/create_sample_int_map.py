import pickle
from sys import argv

def create_sample_int_map(infile):
    sim = {}
    with open(infile, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines[1:]): # ignore header
        name = line.split('\t')[0]
        sim[name] = i
        sim[i] = name
    with open('data/intermediate/sample_int_map.pickle', 'wb') as f:
        pickle.dump(sim, f)

if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError(
            'Usage: python3 create_sample_int_map.py <samples>.tsv')
    infile = argv[1]
    create_sample_int_map(infile)

