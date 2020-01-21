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

1. The first stage in the process is converting a directory full of FASTQ contig
files into a python nested dictionary, and serializing it with pickle. I
already did this. The associated program is pickle_input.py.

1. The second stage is generating a list of unique kmers in the dataset. I used
DSK for that.

1. The third stage is pruning that list down to kmers that appear at least
10 times. The associated program is reduce_input.py.

1. The fourth stage is assembling a database that maps kmers to samples in which
those kmers appear. This allows us to have a contains(sample, kmer) lookup
table. The associated program is create_kmer_db.py. This is extremely
memory-intensive and you can't use a laptop to do it unless you have the power
of resurrection so you can see if it's finished after you die of boredom and
old age a few times. Time scales sub-linearly with available memory because of
the glory of O(1) hash table lookups.

1. The fifth stage is creating a reduced version of that database for PSL. We need
the full thing for the post-psl analysis at the end, but PSL does not need
to know where in the genome a kmer occurs. The associated program is
clean_kmer_db.py.

1. The sixth stage is using the existing data to generate input files for PSL.
The assoicated programs and input files are
    - `convert_kmer_db_to_psl_input.py`, `contains_kmer_sample.txt`
    - `convert_phenos_to_psl_input.py`, `resistance_sample_class.txt`
    - `create_resistance_kmer_class_map.py`, `resistance_kmer_class.txt`
    - `create_similar_antibiotic_map.py`, `similar_antibiotic_antibiotic.txt`
    - `create_similar_sample_map.py`, `similar_sample_sample.txt`

1. The seventh stage is running PSL.

1. The eighth stage is examining the subset of kmers that PSL flags,
and associating them with specific genes, indels, and SNPs.

1. The ninth stage is comparing those attributes to existing datasets,
to determine the accuracy of our method and to hopefully find new genes,
indels, and SNPs that contribute to resistance.

## Installation
Requirements:
    - python 3.6 or later
    - pandas

## Usage
Eventually the entire pipeline will be automated with a script.

## Notes
The full schema for this repository, were you going to run the entire pipeline,
would look something like this. Note that I assume if you have a local DSK
binary, it will be somewhere else.

```
psl-gwas/
    README.md
    contigs/..........FASTQ files containing raw contigs
    data/
        preprocess/...interim files used in different pipeline stages
        psl/..........input files for psl
    src/..............python scripts used in different pipeline stages
    psl/..............psl files defining predicates and input file paths
```

relational psl rule:
    if a kmer is diferent (hamming) from another kmer by exactly 1 base,
    it is likely a SNP/indel (possibility of random also, how to account for that?)

    `KmerResistance(k1) && HammingDiff1(k1, k2) >> !KmerResistance(k2)`

    there could be up to 30 random kmers totally unrelated that have a hamming distance of 1.
    idea is you can identify the pop default, and the mutation holding kmer.

iid method
    - for each kmer, get all kmers that are a hamming distance of 1 away.
       if that kmer is only associated with resistance and another 1d away kmer is not associated with
       resistnace, keep it
    - eliminate all kmers that are present in both resistant and nonresistant samples
    - eliminstae all kmers that are present in nonresistance samples
    score kmers by appearence in diverse samples
    

## Todo
reorganize repo structure
predict_pheno has its own directory containing bin, src, data, psl
gwas has its own directroy containing bin, src, data, psl
preprocess has its own directory containing bin, src, data
    preprocess binaries include scripts for running dsm/dsk

pipeline is now preprocess up to point of gwas and predict pheno divergence
then either gwas or predict pheno preprocess
then either gwas or predict pheno binaries for running psl/iid
