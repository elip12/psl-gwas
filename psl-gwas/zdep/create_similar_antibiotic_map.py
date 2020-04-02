import pandas as pd
import numpy as np
import pickle
from sys import argv

# class int map
def read_cim():
    with open('data/intermediate/class_int_map.pickle', 'rb') as f:
        cim = pickle.load(f)
    return cim

def main(infile):
    df = pd.read_csv(infile, delimiter='\t')
    cim = read_cim()
    
    idcol = df.columns[0]
    df.drop(idcol, axis=1, inplace=True)  
    df.rename(columns=cim, inplace=True) 
    df = df.corr()
    df = df.stack()
    df = df.reset_index()
    # hacky remapping. turn everything under .75 into 0, and remap everything remaining to [0.5,1]
    df[0] = df[0].apply(lambda x: 0 if x < 0.75 else 1 - ((1 - x) * 2))
    df.to_csv('data/preprocessed/similar_antibiotic_antibiotic.txt', sep='\t',
        index=False, header=False)

if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError(
            'Usage: python3 create_similar_antibiotic_map.py <phenotypes>.tsv')
    infile = argv[1]
    main(infile)
