from large_file_processor import main, write_list, parse_args
import pandas as pd
import numpy as np
#from scipy.stats import fisher_exact

def process(data, outfile, THRESH, dfres, dfvul):
    kmer_db_chunk = []
    nphenos = dfres.shape[1]
    for line in data:
        linelist = line.split()
        kmer_db_line = [linelist[0]]
        res = np.zeros(nphenos)
        vul = np.zeros(nphenos)

        for sci in linelist[1:]:
            sample_id = sci.split(',')[0]
            if sample_id not in dfres.index:
                continue
            kmer_db_line.append(sample_id)
            # collect resistant/vulnerable frequencies for each antibiotic for
            # this kmer
            res += dfres.loc[sample_id].to_numpy()
            vul += dfvul.loc[sample_id].to_numpy()

        # 1 test per antibiotic; kmer needs to pass only 1 to avoid
        # getting filtered out
        a = np.where((res + vul >= THRESH) \
                    & (res > 0) \
                    & (vul / res < 0.05))[0]

        if a.size > 0:
            kmer_db_chunk.append(' '.join(kmer_db_line))
    write_list(kmer_db_chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 20
    INPUT_FILE = 'data/intermediate/unitig_sample_map.txt'
    outfile = 'data/intermediate/unitig_sample_map_reduced.txt'
    THRESH = 5 # this value should depend on the min frequency in the phenotype col.
    df = pd.read_csv('data/raw/phenotypes.tsv', delimiter='\t')
    idcol = df.columns[0]
    df.set_index(idcol, inplace=True)
    dfres = df.fillna(0)
    dfvul = 1 - df.fillna(1)
    main(process, NUM_WORKERS, INPUT_FILE, outfile=outfile, THRESH=THRESH, dfres=dfres, dfvul=dfvul)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

