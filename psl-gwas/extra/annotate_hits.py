import pandas as pd
from sys import argv
'''
takes in a reference file (gff3) with start and end locations of each feature,
and an annotations file (sam) with the location of each kmer in the reference
(generated with bwa from a fasta file of the reference),
and cross references them to return the feature each kmer appears in if any.
Likely, most of these will be genes.

Returns kmer, feature pairs with some extra information
'''

def get_name(row):
    atts = row['attributes'].split(';')
    feature_name = 'NA'
    for att in atts:
        if att.startswith('Name='):
            feature_name = att[5:]
    return feature_name

# read reference gff3 file
refnames = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
ref = pd.read_csv(argv[1], sep='\t', skiprows=2, header=None, names=refnames)
ref = ref[ref['feature'] == 'gene']

# read sam annotations file
samnames = ['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual']
sam = pd.read_csv(argv[2], sep='\t', skiprows=2, header=None, names=samnames, usecols=list(range(len(samnames))))

# drop all kmers that are not present in the reference
sam = sam[sam['flag'] != 4]
sam = sam[['qname', 'pos', 'seq']]

# iterate thru kmers and match them to features
n_kmers = int(argv[3])
for i in range(min(n_kmers, sam.shape[0])):
    qname, pos, seq = sam.iloc[i]
    # check if kmer falls within a feature
    within_feature = ref[(ref['start'] <= pos) & (ref['end'] >= pos)]
    for i, row in within_feature.iterrows():
        feature_type = row['feature']
        feature_name = get_name(row)
        print('Query:', qname, 'Seq:', seq, 'Type:', feature_type, 'Name:', feature_name) 
    if len(within_feature) == 0:
        # check if kmer falls between features
        gen = ref.iterrows()
        curr = next(gen)[1]
        prev = None
        while curr['start'] < pos:
            prev = curr
            curr = next(gen)[1]
        prev_feature_type = prev['feature']
        prev_feature_name = get_name(prev)
        curr_feature_type = curr['feature']
        curr_feature_name = get_name(curr)
        print('Query:', qname, 'Seq:', seq, 'Prev Type:', prev_feature_type,
            'Prev Name:', prev_feature_name, 'Next Type:', curr_feature_type,
            'Next Name:', curr_feature_name)

