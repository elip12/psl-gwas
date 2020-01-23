import numpy as np
from sklearn.svm import LinearSVC
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix
from random import shuffle

'''
convert kmer sample dict into matrix with samples as rows, kmers as columns.
presence or absence of kmer is value for each cell
'''

dpath = 'preprocess/data/preprocessed'

size = 15000000
X = np.zeros((300, size), dtype=np.ushort)
y = np.zeros(300)
test_X = np.zeros((55, size), dtype=np.ushort)
test_y = np.zeros(55)

samples = list(range(355))
shuffle(samples)
train = samples[:300]
heldout = samples[300:]

# fill matrix
# lines is list of all lines in file
with open(f'{dpath}/contains_sample_kmer.txt', 'r') as f:
    lines = f.readlines()
for l in lines:
    r,c,v = l.split('\t')
    r = int(r)
    c = int(c)
    v = np.ushort(int(float(v)))
    if r in heldout:
        test_X[heldout.index(r)][c] = v
    else:
        X[train.index(r)][c] = v

# lines is now a list of resistances
with open(f'{dpath}/resistance_sample_class.txt', 'r') as f:
    lines = f.readlines()
for l in lines[332:661]:
    r, _, v = l.split('\t')
    r = int(r)
    v = np.ushort(int(float(v)))
    if r in heldout:
        test_y[heldout.index(r)] = v
    else:
        y[train.index(r)] = v


svm = LinearSVC()
svm.fit(X,y)
predicted = svm.predict(test_X)
print(accuracy_score(test_y, predicted))
print(precision_recall_fscore_support(test_y, predicted))
print(confusion_matrix(test_y, predicted))
print(predicted)
print()
print(test_y)
