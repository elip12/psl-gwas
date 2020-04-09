from sys import argv

gene = argv[1]
abfile = argv[2]

with open(abfile, 'r') as f:
    lines = f.readlines()
seq = ''
found = False
for line in lines:
    if line.startswith(f'>{gene}'):
        found = True
        continue
    if found == True and line.startswith(f'>'):
        break
    if found == True:
        seq += line.rstrip()
with open('seq.fa', 'w') as f:
    f.write(seq)

