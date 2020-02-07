from large_file_processor import main, write_list, parse_args
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

def process(data, outfile, THRESH, df):
    kmer_db_chunk = []
    for line in data:
        linelist = line.split()
        kmer_db_line = [linelist[0]]
        prop_res = [[0,0] for _ in range(df.shape[1])]
        for sci in linelist[1:]:
            sample_id = sci.split(',')[0]
            if sample_id not in df.index:
                continue
            kmer_db_line.append(sample_id)
            # collect resistant/vulnerable frequencies for each antibiotic for
            # this kmer
            slice_ = df.loc[sample_id]
            for i, antibiotic in enumerate(df.columns):
                val = slice_[antibiotic]
                if val == 1.0:
                    prop_res[i][0] += 1
                elif val == 0.0:
                    prop_res[i][1] += 1
                # if val is NA, do nothing
        # 1 chi-squared test per antibiotic; kmer needs to pass only 1 to avoid
        # getting filtered out
        kmer_passed = False
        for observed in prop_res:
            if sum(observed) < THRESH:
                continue
            # null hypothesis in this test is that this kmer connotes resistance
            #p = fisher_exact(table)[1]
            
            # this is kinda iffy but we are doing it for now
            p = 1 if observed[0] == 0 else observed[1] / observed[0]
            if p < 0.05:
                # print(p, observed, linelist[0])
                kmer_passed = True
                break
        if kmer_passed:
            kmer_db_chunk.append(' '.join(kmer_db_line))
    write_list(kmer_db_chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 1
    INPUT_FILE = 'data/intermediate/kmer_sample_map_head.txt'
    outfile = 'data/intermediate/kmer_sample_map_reduced_test.txt'
    THRESH = 5 # this value should depend on the min frequency in the phenotype col.
    df = pd.read_csv('data/intermediate/abr_resist_phenos.tsv', delimiter='\t')
    df.drop(['Date', 'Species', 'Tissue'], axis=1, inplace=True)
    df.set_index('Sample', inplace=True)
    main(process, NUM_WORKERS, INPUT_FILE, outfile=outfile, THRESH=THRESH, df=df)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

