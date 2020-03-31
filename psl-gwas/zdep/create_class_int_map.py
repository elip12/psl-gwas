import pickle
import pandas as pd
from sys import argv

def create_class_int_map(fname):
    df = pd.read_csv(fname, sep='\t')
    idcol = df.columns[0]
    df.drop(idcol, axis=1, inplace=True)
    cim = {}
    for i, c in enumerate(df.columns.values):
        cim[i] = c
        cim[c] = i
    with open('data/intermediate/class_int_map.pickle', 'wb') as f:
        pickle.dump(cim, f)

if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError(
            'Usage: python3 create_class_int_map.py <phenotypes>.tsv')
    infile = argv[1]
    create_class_int_map(infile)

