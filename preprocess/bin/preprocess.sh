#clean_kmer_db.py                 create_kmer_db.py                    create_sample_int_map.py          __init__.py              __pycache__
#convert_kmer_db_to_psl_input.py  create_kmer_int_map.py               create_similar_antibiotic_map.py  large_file_processor.py  reduce_input.py
#convert_phenos_to_psl_input.py   create_resistance_kmer_class_map.py  create_similar_sample_map.py      pickle_input.py

python3 src/reduce_input.py
python3 src/create_kmer_db.py
python3 src/clean_kmer_db.py
python3 src/create_resistance_kmer_class_map.py
python3 src/create_sample_int_map.py
python3 src/create_kmer_int_map.py
python3 src/create_similar_antibiotic_map.py
python3 src/create_similar_sample_map.py
python3 src/convert_kmer_db_to_psl_input.py
python3 src/convert_phenos_to_psl_input.py
