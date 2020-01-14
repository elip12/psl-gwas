import numpy as np
from sklearn.svm import LinearSVC
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix

'''
convert kmer sample dict into matrix with samples as rows, kmers as columns.
presence or absence of kmer is value for each cell
'''

X = np.zeros((300, 9176498))
y = np.zeros(300)
test_X = np.zeros((55, 9176498))
test_y = np.zeros(55)

# fill matrix
# lines is list of all lines in file
with open('../../data/psl/contains_sample_kmer_head.txt', 'r') as f:
    lines = f.readlines()
for l in lines:
    r,c,v = l.split('\t')
    if int(r) >= 300:
        test_X[int(r) - 300][int(c)] = float(v)
    else:
        X[int(r)][int(c)] = float(v)

# lines is now a list of resistances
with open('../../data/psl/resistance_sample_class.txt', 'r') as f:
    lines = f.readlines()
for l in lines[:331]: # this will change
    r, _, v = l.split('\t')
    if int(r) >= 300:
        test_y[int(r) - 300] = float(v)
    else:
        y[int(r)] = float(v)


svm = LinearSVC()
svm.fit(X,y)
predicted = svm.predict(test_X)
a = accuracy_score(test_y, predicted)
print(a)
