infile = 'data/postprocessed/scored_kmers.txt'
outfile = 'data/postprocessed/scored_kmers.fsa'

with open(infile, 'r') as f:
    lines = f.readlines()
lines = [line.split('\t')[0] for line in lines]
out = [f'>{i}\n{line}' for i, line in enumerate(lines)]

with open(outfile, 'w') as f:
    f.writelines('\n'.join(out))
