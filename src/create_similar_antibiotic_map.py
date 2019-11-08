import pandas as pd
import numpy as np
import pickle

# class int map
def read_cim():
    with open('data/preprocess/class_int_map.pickle', 'rb') as f:
        cim = pickle.load(f)
    return cim

def main():
    df = pd.read_csv('data/preprocess/abr_resist_phenos.tsv', delimiter='\t')
    cim = read_cim()
    
    df.drop(['Sample', 'Date', 'Species', 'Tissue'], axis=1, inplace=True)  
    df.rename(columns=cim, inplace=True) 
    df = df.corr()
    df = df.stack()
    df = df.reset_index()
    df[0] = min(df[0], 0)
    df.to_csv('data/psl/similar_antibiotic_antibiotic.txt', sep='\t',
        index=False, header=False)

if __name__ == '__main__':
    main()

