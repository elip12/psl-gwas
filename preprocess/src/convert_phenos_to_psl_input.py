import pandas as pd
import pickle
from sys import argv

ipath = 'data/intermediate'
ppath = 'data/preprocessed'

# sample int map
def read_sim():
    with open(f'{ipath}/sample_int_map.pickle', 'rb') as f:
        sim = pickle.load(f)
    return sim

def read_cim():
    with open(f'{ipath}/class_int_map.pickle', 'rb') as f:
        cim = pickle.load(f)
    return cim

def main(infile):
    df = pd.read_csv(infile, delimiter='\t')
    sim = read_sim()
    cim = read_cim()
    
    idcol = df.columns[0]
    df['sample_id'] = df[idcol].apply(lambda x: sim[x])
    df.drop(idcol, axis=1, inplace=True)
    df.rename(columns=cim, inplace=True) 
    
    new = pd.DataFrame(columns=['sample_id', 'class', 'value'])
    for col in df:
        if col != 'sample_id':
            slice_ = df[['sample_id', col]].copy()
            slice_['class'] = col
            slice_.rename(columns={col: 'value'}, inplace=True)
            new = pd.concat([new, slice_[['sample_id', 'value', 'class']]],
                ignore_index=True, sort=False) 
    df = new.dropna()
    df.to_csv(f'{ppath}/resistance_sample_class.txt', sep='\t',
        index=False, header=False)

if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError(
            'Usage: python3 convert_phenos_to_psl_input.py phenotypes.tsv')
    infile = argv[1]
    main(infile)

