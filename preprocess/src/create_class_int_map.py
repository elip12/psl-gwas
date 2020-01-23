import pickle
import pandas as pd

def create_class_int_map():
    kim = {}
    fname = '/data/intermediate/abr_resist_phenos.tsv'
    df = pd.read_csv(fname)
    df.drop(['Sample', 'Date', 'Species', 'Tissue'], axis=1, inplace=True)
    cim = {}
    for i, c in enumerate(df.columns.values):
        cim[i] = c
        cim[c] = i
    with open('data/intermediate/class_int_map.pickle', 'wb+') as f:
        pickle.dump(cim, f)

if __name__ == '__main__':
    create_class_int_map()

