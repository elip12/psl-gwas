


base = 'ecoli/data/preprocessed/'
truths_file = base + 'truth_unitig_pheno.txt'
unitig_int_map = load_pickle(base + 'unitig_int_map.pkl')
unitig_sample_map_file = base + 'unitig_sample_map.txt'
phenos_df = pd.read_csv('ecoli/data/raw/phenos.tsv', sep='\t')
sim = load_pickle(base + 'sample_int_map.pkl')
pim = load_pickle(base + 'pheno_int_map.pkl')

def get_truth_kmer_seqs_that_occur_in_data():
    with open(truths_file, 'r') as f:
        truthlines = f.readlines()
    #with open(unique_kmers_file, 'r') as f:
    #    kmerslines = f.readlines()
    #kmerslines = set(k.split('\t')[0] for k in kmerslines)
    truthslines = [(unitig_int_map[int(k)], p) for k,p,_ in truthslines.split('\t')]
    return truthslines

def read_unitig_sample_map_into_dict():
    with open(unitig_sample_map_file, 'r') as f:
        lines = f.readlines()
    lines = {l.split('\t')[0]:l.split('\t')[1:] for l in lines}
    return lines

def cross_truth_seqs_with_usm(usm, truthslines):
    return [(k,p, usm[k]) for k, p in truthslines]

def convert_sample_ids_to_pheno_value(data)
    means = []
    for k, p, samples in data:
        samplelist = [phenos_df[pim[p]][sim[s]] for s in samples]
        samplelist = [s for s in samplelist if s != 'NA']
        mean = np.mean(samplelist)
        means.append(mean)
    print('min: ', min(means))
    print('q1:  ', np.percentile(means, 25))
    print('med: ', np.percentile(means, 50))
    print('mean:', np.mean(means))
    print('q3:  ', np.percentile(means, 75))
    print('max: ', max(means))

def main():
    truths = get_truth_kmer_seqs_that_occur_in_data()
    usm = read_unitig_sample_map_into_dict()
    crossed = cross_truth_seqs_with_usm(usm, truths)
    convert_sample_ids_to_pheno_values(crossed)

