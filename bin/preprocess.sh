#!/bin/bash
python3 psl-gwas/count_kmers.py "$@" \
&& python3 psl-gwas/preprocess.py "$@" \
&& python3 psl-gwas/psl_input.py "$@"

