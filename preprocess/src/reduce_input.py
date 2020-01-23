# removes all kmers that occur fewer than 10 times in the dataset.
def reduce_input():
with open('data/intermediate/unique_kmers.txt', 'r') as fin:
    with open('data/intermediate/unique_kmers_reduced.txt', 'w+') as fout:
        for line in fin:
            s, n = line.split()
            if int(n) >= 10:
                fout.write(line)

if __name__ == '__main__':
    reduce_input()

