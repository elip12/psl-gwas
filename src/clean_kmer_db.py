from large_file_processor import main, write_list, parse_args

def process(data, outfile, THRESH, df):
    kmer_db_chunk = []
    for line in data:
        linelist = line.split()
        if len(linelist) <= THRESH:
            continue
        kmer_db_line = [linelist[0]]
        prop_res = [0 for _ in range(df.shape[1])]
        for sci in linelist[1:]:
            sample_id = sci.split(',')[0]
            kmer_db_line.append(sample_id)
            # replace this part with a X-squared test maybe
            for i, antibiotic in enumerate(df):
                if df[antibiotic][sample_id] == '1':
                    prop_res[i] += 1
        prop_res = [nr / (len(linelist) - 1) for nr in prop_resistant]
        kmer_passed = False
        for pr in prop_res:
            # if more than 5% of samples that contain this kmer are not resistant,
            # this kmer does not confer resistance. The 5% is pretty arbitrary,
            # and is for the possibility of epistatic interactions
            if pr >= 0.95:
                kmer_passed = True
                break
        if kmer_passed:
            kmer_db_chunk.append(' '.join(kmer_db_line))
    write_list(kmer_db_chunk, outfile)

def main_wrapper():
    NUM_WORKERS = 16
    INPUT_FILE = 'data/preprocess/kmer_sample_map.txt'
    outfile = 'data/preprocess/kmer_sample_map_reduced.txt'
    THRESH = 10
    df = pd.read_csv('data/preprocess/abr_resist_phenos.tsv', delimiter='\t') 
    df.set_index('Sample', inplace=True)
    main(process, NUM_WORKERS, INPUT_FILE, outfile=outfile, THRESH=THRESH, df=df)

if __name__ == '__main__':
    parse_args()
    main_wrapper()

