# PSL GWAS
Eli Pandolfo

## Description
This project aims to perform a genome-wide association study (GWAS) on microbial
genomes (specifically E. coli), using Probabalistic Soft Logic (PSL).

A GWAS is a computational method of mapping specific attributes of a genotype
(such as SNPs, indels, or genes) to specific attributes of a phenotype (in this
case, resistance to different types of antibiotics). Formulated as a machine
learning problem, our model uses Bayesian inference to predict which genotype
attributes are associated with resistance to certain antibiotics, given
the presence or absence of that genotype attribute in a sample that either has
or lacks resistance to those antibiotics. We use kmers of length 30 to represent
genotype attributes, since they capture both point mutations and entire genes,
and are usable with unaligned genomes (we have only contigs of length 1000 to
100000 bps).

## Repo structure
-   *predict-pheno*
    Train and test a logistic regression classifier to predict the phenotype of a
    sample.
-   *preprocess*
    Python and bash scripts for the preprocessing pipeline
-   *postprocess*
    Python and bash scripts for the postprocessing pipeline

## Installation
1. create virtual environment for python modules and activate it
1. `./preprocess/bin/setup.sh`
This will install all necessary python packages, and create the `data/`
directory and subdirectories, which will not be tracked by git.

## Usage
1.  Create phenotypes file. This should be a tab-separated file of the form
    `ID    pheno1    pheno2...`
    where ID is the unique string sample ID (unique among other samples)
    and pheno* is in {1, 0, NA}, meaning a sample displays a phenotype,
    does not display a phenotype, or pheno information is not known.
    Put this file, which should be end in `.tsv`, in `data/raw/`
1.  Create samples file. This should be a tab-separated file of the form
    `ID    path`
    where ID is the unique string sample ID (corresponding to IDs in
    the phenotype file) and path is the absolute or relative path from the root
    of the repository to the `.fa` file holding reads for that sample.
    Put this file, which should end in `.tsv`, in `data/raw`.
    It is fine for the FASTA files themselves to be in `data/raw/`,
    but not necessary.
1.  more to follow...



