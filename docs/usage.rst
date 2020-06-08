Usage
#####

Quick Start
===========

``./bin/run.sh --project <project> --sample <samples.tsv> --pheno <phenos.tsv>``

Detailed Usage
==============

The global run script (``bin/run.sh``) is idempotent, meaning the program fails,
you can rerun it with the same arguments and it will skip any preprocessing
files that already exist. Note that it is possible for the program to fail
having written some, but not all data to a file. In this case, you will need
to delete the file (in ``<project>/data/preprocessed``) so that PSL-GWAS
does not skip that file.

PSL-GWAS automatically uses the max number of available threads and 95% of the max
memory of the system it runs on.

PSL-GWAS takes in *de novo* assembled contigs or fully sequenced genomes
in FASTA (.fa, .fsa) format.

#. Create a file holding the string IDs and paths to your sample genomes.
    This is a tab-separated file with two columns. It must include a header but
    the column names in the header can be anything. The first column holds
    unique string IDs for each sample, and the second column holds paths
    to that sample's FASTA file. The paths are should either be absolute or
    relative to the root of the repo. So, if you store your sample genomes
    in a directory called ``contigs/`` which you put in the ``raw/`` directory,
    your sample file would look like this:

    =====   =====  
    ID      Path      
    =====   =====  
    s1      <project>/data/raw/contigs/s1_ctgs.fa
    s2      <project>/data/raw/contigs/s2_ctgs.fa 
    s3      <project>/data/raw/contigs/s3_ctgs.fa
    ...     ...
    =====   =====

    You can call this file anything. ``samples.tsv`` is simple and easy.
    This file goes in ``<project>/data/raw/``


#. Create a file holding the string IDs and phenotype values for your samples.
    This is a tab-separated file with at least two. It must include a header.
    The first column of the header can be anything (we recommend "ID"), and the
    remaining columns of the header should be the names of the phenotypes
    you are testing.
    The first column holds unique string IDs for each sample (that must
    match the ids in the samples file), and the remaining columns hold
    data on whether each sample displays the ID. When phenotype information
    is not known, use "NA". Phenotype information can be binary {0, 1}
    or continuous [0, 1]. 0 means the sample does not display the phenotype,
    1 means it displays it with high strength, and intermediate values mean
    it displays it with weaker strength.

    If you are testing antibiotic resistance, your phenotypes file could look
    like this:

    =====   ==========  ==========  ===
    ID      Penicillin  Ampicillin  ...
    =====   ==========  ==========  ===
    s1      0           0.75        ...
    s2      1           NA          ...
    s3      1           1           ...
    ...     ...         ...         ...
    =====   ==========  ==========  ===

    You can call this file anything. ``phenos.tsv`` is simple and easy.
    This file goes in ``<project>/data/raw/``

#. Run the gwas from the root of the repo
    ``./bin/run.sh --project <project> --sample <samples.tsv> --pheno <phenos.tsv> ...``

Options Reference for run.sh
============================

Required flags:
    --project PROJECT     Name of project, defined with startproject.sh <project>.
    --sample SAMPLE       Basename of samples file. Ex: ``samples.tsv``
    --pheno PHENO         Basename of phenos file. Ex: ``phenos.tsv``
      
Optional flags:
    -d, --debug           More verbose logging.
    -k K, --k K           Kmer length in nucleotide bases. Default: 31
    --minkf MINKF         Minimum kmer frequency used during preprocessing filtering. Default: 0.01
    --maxkf MAXKF         Maximum kmer frequency, used during preprocessing filtering. Default: 0.95
    --correlation-thresh CORRELATION_THRESH
                          Penetrance threshold for filtering in preprocessing. Default: 0.5
    -p, --param           Ignore k, minkf, maxkf, and correlation-thresh options
                            and use param file in project directory.
    --truth TRUTH         Fasta file holding truths data for benchmarking or weight learning.
                          Labels correspond to phenos, sequences hold genes or
                          unitigs that cause the phenotype.
    --postgres DATABASE   Postgres database the default user (usually your username) can access.    

