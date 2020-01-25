import random

#this will split data/resistance_sample_class_ab1.txt into 
#...ab1_train.txt, ...ab1_truth.txt, and create an partial ...ab1_target.txt.

inpath = 'data/preprocessed'
outpath = 'predict-pheno/data'
with open(f'{inpath}/resistance_sample_class_ab1.txt', 'r') as f:
    lines = f.readlines()
random.shuffle(lines)
length = len(lines)
cutoff = int(0.8 * length)
reduce_ = lambda x: '\t'.join((x.split('\t')[0], x.split('\t')[2]))
train_data = [reduce_(line) for i, line in enumerate(lines) if i < cutoff]
test_data = [reduce_(line) for i, line in enumerate(lines) if i >= cutoff]
with open(f'{outpath}/resistance_sample_class_ab1_train.txt', 'w') as f:
    f.writelines(train_data)
with open(f'{outpath}/resistance_sample_class_ab1_test.txt', 'w') as f:
    f.writelines(test_data)
        
