from sys import argv

fname = argv[1]
outfile = 'most_common_ab_genes.out'
with open(fname, 'r') as f:
    lines = f.readlines()
maxbits = 0
maxgene = None
for line in lines:
    if line.startswith('Sequence') or line.startswith('\n') or line.startswith('--'):
        continue
    gene, bits, _ = line.rstrip().split()
    if float(bits) > maxbits:
        maxbits = float(bits)
        maxgene = gene
with open(outfile, 'a+') as f:
    f.write(maxgene + '\n')
