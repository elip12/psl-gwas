import pandas as pd
from sys import argv

df = pd.read_csv(argv[1], sep='\t', header=None)
for col in df[1].unique():
    dfcol = df[df[1] == col]
    dfcol = dfcol.nlargest(int(argv[3]), columns=2)
    if col == 'Gentamycin':
        col = 'Gentamicin'
    dfcol.to_csv(f'{argv[2]}/psl.{col.lower()}', sep='\t')
