with open('pyseer.assoc', 'r') as f:
    lines = f.readlines()
lines = [line.split('\t')[0] for line in lines]
out = [f'>{i}\n{line}' for i, line in enumerate(lines)]
with open('pyseer_assoc.fsa', 'w') as f:
    f.writelines('\n'.join(out))
