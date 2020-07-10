#!/usr/bin/env python3 
from sys import argv

with open(argv[1], 'r') as f:
    lines = f.readlines()
with open(argv[2], 'w') as f:
    for i, line in enumerate(lines):
        f.write(f'>{i}\n{line}')

