INDIR='/Users/Eli/Documents/school/spring20/linqs/raw/contigs'
ABDB='/Users/Eli/Documents/school/spring20/linqs/blast/resfinder_db/abclasses.fsa'
SDB='/Users/Eli/Documents/school/spring20/linqs/blast/samplesdb/samples.db'
find_most_common() {
    for fname in $INDIR/*.fa; do
        f="${fname##*/}"
        blastn -db "$ABDB" -query "$INDIR/$f" -out "f.out"
        grep -A 2 "Sequences producing" f.out > f_consolidated.out
        python3 get_highest_bits.py f_consolidated.out
    done
    rm f.out
    rm f_consolidated.out
    sort most_common_ab_genes.out | uniq -c | sort -r  > genes.out
    rm most_common_ab_genes.out
}

get_sample() {
    cat genes.out | while read line; do
        gene=$(echo $line | cut -d ' ' -f 2)
        python3 get_gene_seq.py $gene $ABDB
        echo ">$gene" >> genes.fa
        cat "seq.fa" >> genes.fa
        echo >> genes.fa
        blastn -db "$SDB" -query "seq.fa" -out "seq.out"
        python3 get_all_samples.py seq.out $gene
    done
    rm seq.fa
    rm seq.out
}

create_phenos() {
    python3 create_phenos_file.py
    rm samples_containing_gene_*
    rm genes.out
}

main() {
    find_most_common
    get_sample
    create_phenos
}

main

