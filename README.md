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

Contact Eli Pandolfo (epandolf@ucsc.edu) with questions.

## Repo structure
`bin`: holds global shell scripts for initializing and running gwas
`psl-gwas`: holds python and psl source files
`README.md`: this file
`requirements.txt`: python requirements

## Installation
1. create virtual environment for python modules and activate it
1. `cd` into your clone of this repo
1. `./bin/startproject.sh <myproject>`
This will start a psl-gwas project, and create a directory for your project.
It will create the necessary files and directories in your project directory.

## Usage
1.  Create phenotypes file. This should be a tab-separated file of the form
    `ID    <pheno1>    <pheno2>...`
    where ID is the unique string sample ID (unique among other samples)
    and the values for each phenotype is in {NA, 0, int}, meaning pheno info
    is not known, a sample does not display a phenotype, or a sample displays
    a phenotype to some varying degree. You should use the real names of
    the phenotypes in this file. It does not matter what the first column
    is explicitly called, but the first column needs to hold the IDs.
    Put this file, which should be end in `.tsv`, in `<myproject>/data/raw/`
1.  Create samples file. This should be a tab-separated file of the form
    `ID    path`
    where ID is the unique string sample ID (corresponding to IDs in
    the phenotype file) and path is the absolute or relative path from the root
    of the repository to the `.fa/.fsa` file holding reads for that sample.
    Put this file, which should end in `.tsv`, in `<myproject>/data/raw`.
    It is fine for the FASTA files themselves to be in `<myproject>/data/raw/`,
    but not necessary.
1.  Run the gwas:
    `./bin/run.sh --project <myproject> --samples <samples> --phenos <phenos> ...
    TODO: more guidance on running.
    The processing pipeline is fully idempotent, so if you cannot run the full
    pipeline all the way through, you can do `./bin/run.sh ...` again and it
    will start where you left off.

After it runs: in `<myproject>/data/postprocessed/`, there will be a file called
`scored_unitigs.fsa`. This is a fasta file holding the unitigs with > 95% confidence
to any phenotype.

## Acknowledgements
Thanks to the UCSC LINQS Lab and the UCSC Camps lab.

## License
This project is licensed under the MIT license.

## Future
Ability to download test data and replicate results.

