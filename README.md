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

