from sys import argv
fname = argv[1]
gene = argv[2]
outfile = f'samples_containing_gene_{gene}.out'
all_samples = '../../raw/samples.tsv'

with open(all_samples, 'r') as f:
    lines = f.readlines()
all_samples = [line.split('\t')[0] for line in lines][1:]

with open(fname, 'r') as f:
    lines = f.readlines()
samples = set()

for line in lines:
    for sample in all_samples:
        if line.startswith(sample):
            samples.add(sample + '\n')
with open(outfile, 'a+') as f:
    f.writelines(list(samples))
