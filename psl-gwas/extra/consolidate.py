from sys import argv
from os.path import basename
import pandas as pd

"""
This script takes in a file holding newline separated kmers.
It loops through all pairs of kmers, and if two kmers share a significant
portion of their sequence, it consolidates them into a single kmer.

Max value of K is 30, since no kmer will be shorter than 30 kmers
Min value of K I set to 20, which seems to be a nice number. You can
check how many kmers get consolidated at each value of k, and for
ceftazidime, there is no difference between min k of 22 and min k of 10.
Anything under 10 has a high enough probability of randomly matching that I
dont want to look at it.
"""

def consolidate(data, k, p=False):
    for i, (unitig1, ranks1) in enumerate(data):
        for j in range(i + 1, len(data)):
            unitig2, ranks2 = data[j]
            if unitig1[:k] == unitig2[-k:]:
                newunitig = unitig2[:-k] + unitig1
                if p:
                    print(j, '\t', unitig2)
                    print(i, '\t', unitig1.rjust(len(newunitig), ' '))
                    print('-\t', '-'.rjust(len(newunitig), '-'))
                    print(j, '\t', newunitig, '\n')
                data[i] = None
                data[j] = (newunitig, ranks1 + ranks2)
                break
            elif unitig1[-k:] == unitig2[:k]:
                newunitig = unitig1[:-k] + unitig2
                if p:
                    print(i, '\t', unitig1)
                    print(j, '\t', unitig2.rjust(len(newunitig), ' '))
                    print('-\t', '-'.rjust(len(newunitig), '-'))
                    print(j, '\t', newunitig, '\n')
                data[j] = (newunitig, ranks1 + ranks2)
                data[i] = None
                break
    return [d for d in data if d is not None]


def main():
    name = basename(argv[1])
    unitigs = pd.read_csv(argv[1], sep='\t', index_col=0)
    unitigs = list(unitigs.iloc[:, 0])
    unitigs = [(u, [1/(i+1)]) for i, u in enumerate(unitigs)]
    print('Original num kmers:', len(unitigs))
    for k in range(30, 20, -1):
        unitigs = consolidate(unitigs, k, p=False)
        #print(len(unitigs))
    print('Final num kmers:', len(unitigs))
    unitigs = [(unitig, 1 / (sum(ranks) / len(ranks))) for unitig, ranks in unitigs]
    unitigs = sorted(unitigs, key=lambda x: x[1])
    unitigs = [u[0] for u in unitigs]
    with open(f'{argv[2]}/{name}', 'w') as f:
        for u in unitigs:
            f.write(u + '\n')
        #print(u)


main()
