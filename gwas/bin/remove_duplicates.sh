#perl -i.bak -ne 'print if ! $x{$_}++' data/preprocessed/contains_sample_kmer.txt
perl -i.bak -ne 'print if ! $x{$_}++' data/preprocessed/resistance_kmer_class.txt
