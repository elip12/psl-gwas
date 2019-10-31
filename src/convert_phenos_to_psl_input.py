import pandas as pd
import pickle

# sample int map
def read_sim():
    with open('data/preprocess/sample_int_map.pickle', 'rb') as f:
        sim = pickle.load(f)
    return sim

def read_cim():
    with open('data/preprocess/class_int_map.pickle', 'rb') as f:
        cim = pickle.load(f)
    return cim

def main():
    df = pd.read_csv('data/preprocess/abr_resist_phenos.tsv', delimiter='\t')
    sim = read_sim()
    cim = read_cim()
    
    df.drop(['Sample', 'Date', 'Species', 'Tissue'], axis=1, inplace=True)  
    df['sample_id'] = df['Sample'].apply(lambda x: sim[x])
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
    df.to_csv('data/psl/resistance_sample_class.txt', sep='\t',
        index=False, header=False)

if __name__ == '__main__':
    main()

