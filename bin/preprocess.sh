#!/bin/bash
###############################################################################
##  preprocess.sh
##  Invokes preprocessing py scripts, passing all command-line arguments
##  through. run.sh will parse the arguments and only give this script the
##  ones it needs.
##  -   count_kmers.py creates a tab-separated file called unique_kmers.txt,
##      holding all uniique kmers and their counts, similar to DSK.
##  -   preprocess.py maps each kmer to the samples it appears in, samples
##      a random subset of those kmers to build a sample similarity matrix,
##      consolidates kmers into longer unitigs as much as possible without 
##      losing any information, and filters the unitigs based on whether they
##      appear in samples that do not display the pheno.
##  -   pslprep.py converts the preprocessed data into the form psl takes as
##      input.
###############################################################################
CONTAINS='data/preprocessed/contains_obs.txt'
TARGET='data/preprocessed/unitigPheno_target.txt'
BLOCK='data/preprocessed/block_obs.txt'
python3 psl_gwas/count_kmers.py "$@" \
&& python3 psl_gwas/preprocess.py "$@" \
&& python3 psl_gwas/pslprep.py "$@" \
&& cp "$2/$TARGET" "$2/$BLOCK"

