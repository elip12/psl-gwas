from sys import argv

# removes all kmers that occur fewer than 10 times in the dataset.
def reduce_input(n):
    with open('data/intermediate/unique_kmers.txt', 'r') as fin:
        inlines = fin.readlines()
    upper = int(0.98 * n)
    lower = int(0.02 * n)
    test = lambda x: x <= upper and x >= lower
    outlines = [line for line in inlines if test(int(line.split('\t')[1]))]
    with open('data/intermediate/unique_kmers_reduced.txt', 'w') as fout:
        fout.writelines(outlines)

if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError(
            'Usage: python3 reduce_input.py <samples>.tsv')
    infile = argv[1]
    with open(infile, 'r') as f:
        n = len(f.readlines()) - 1
    reduce_input(n)

