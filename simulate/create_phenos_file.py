import pandas as pd
from os import listdir

all_samples = '../raw/samples.tsv'
with open(all_samples, 'r') as f:
    lines = f.readlines()
all_samples = [line.split('\t')[0] for line in lines][1:]

all_genes = './genes.out'
with open(all_genes, 'r') as f:
    lines = f.readlines()
all_genes = [line.split()[1] for line in lines]

df = pd.DataFrame(columns=['ID', *all_genes])
df.set_index('ID', inplace=True)
zeros = [0 for _ in all_genes]
for sample in all_samples:
    df.loc[sample] = zeros

startsterm = 'samples_containing_gene_'
endsterm = '.out' 
for fname in listdir():
    if fname.startswith(startsterm):
        gene = fname[len(startsterm):-len(endsterm)]
        with open(fname, 'r') as f:
            lines = f.readlines()
        for sample in lines:
            df[gene][sample.rstrip()] = 1
df.to_csv('phenos.tsv', sep='\t')
        
