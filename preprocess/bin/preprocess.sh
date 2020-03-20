# $1 is the samples tsv
# $2 is the phenotypes tsv

./preprocess/bin/cat_samples.sh $1

# these take negligible time compared to other parts of the pipeline
# maps each phenotype to a unique int
python3 preprocess/src/create_class_int_map.py $2

# creates PSL input file for SIMILARANTIBIOTIC
python3 preprocess/src/create_similar_antibiotic_map.py $2

# maps each sample to a unique int
python3 preprocess/src/create_sample_int_map.py $1

# creates PSL input file for SAMPLERESISTANCE
python3 preprocess/src/convert_phenos_to_psl_input.py $2

# these all maximize the resources of the server
# count all kmers in input FASTA files
python3 preprocess/src/count_kmers.py

# remove all kmers that appear too many or two few times
python3 preprocess/src/reduce_input.py $1

# pickle input into python dict for portability
python3 preprocess/src/pickle_input.py $1

# map those kmers back to the samples they appear in
python3 preprocess/src/create_kmer_db.py

# take random sample of kmers to input data for PSL rule SIMILARSAMPLE
python3 preprocess/src/create_similar_sample_map.py

# reduce kmer database with simple association test
python3 preprocess/src/clean_kmer_db.py

# map each remaining kmer to a unique int
python3 preprocess/src/create_kmer_int_map.py

# create input data for PSL rule KMERRESISTANCE (target)
python3 preprocess/src/create_resistance_kmer_class_map.py

# create input data for PSL rule CONTAINS
python3 preprocess/src/convert_kmer_db_to_psl_input.py

# remove all duplicate kmers that end up in CONTAINS and KMERRESISTANCE
./preprocess/bin/remove_duplicates.sh
