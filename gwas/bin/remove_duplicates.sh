perl -i.bak -ne 'print if ! $x{$_}++' data/processed/contains_sample_kmer.txt
