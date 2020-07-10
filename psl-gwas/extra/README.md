# annotate_fna_with_sample.py
takes in a bunch of fna files from prodigal output, and combines them into a
single long fasta file so that you can blast your kmers against it then use parse
blast output to get the annotaions for each kmer.
# convert_to_fa.py
converts the output of psl-gwas into fasta format, where annotaitons for each
kmer are just the rank (starting at 0)
# eval_truths.py
when you know the truth data, gives metrics on how many truth kmers were
thrown out in preprocessing. used to tune correlation, minkf, maxkf values.
# parse_blast_output.py
when you blast psl-gwas output kmers against some blast database, this script
returns a list of all hits for each kmer
# annotate_hits.py
takes in a gff3 file with a refernce sequence, and a sam file created by
annotating a fasta-converted kmers output from psl gwas with a fasta format
of the reference with bwa. combines them to get the gene/feature each kmer is
associated with in the reference. TODO: consolidate pipeline.
# disperse.sh
deprecated
# evaluate.py
given some truth data (in fasta format, where annotations hold the name of the
truth sequence and then underscore-separated associated phenotypes), and a
pyseer output holding ranked kmers (only), determines the reciprocal rank,
precision, and recall for that file. used in conjunction with disperse to get
different files with the top 10, 50, etc. TODO: integrate num kmers into this.
# pyseer_to_baseline.py
converts the output of pyseer to a baseline fa format where the annotations
are formatted specifically to be interpreted by pyseer and hold the phenotype
and the lrt-pvalue converted to a confidence score
