from utility import process_file, write_dict, parse_args, printd, \
check_outfile, get_params, write_list
from multiprocessing import Queue, Manager
from collections import Counter
from os import remove
from os.path import join

# cat all samples together
def cat_samples(samples, outfile):
    output_lines = []
    with open(samples, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        sample, fname = line.rstrip().split('\t')
        if not (fname.endswith('.fa') or fname.endswith('.fsa')):
            continue
        with open(fname, 'r') as f:
            write_list(f.readlines(), outfile)
   return seqs 

# takes in a chunk of contigs, and creates a counter holding all the kmers
# and their counts in that chunk
def process(data, q, k):
    counter = Counter()
    for line in data:
        if line.startswith('>') or len(line) < k:
            continue
        for i in range(len(line) - k):
            kmer = line[i: i + k]
            counter[kmer] += 1
    q.put(counter)

def main():
    # load params 
    params = get_params()
    project = params['project']
    k = params['k']

    samples_file = join(project, 'data', 'raw', params['samples'])
    outfile = join(project, 'data', 'preprocessed', 'unique_kmers.txt')
    catted_samples = join(project, 'data', 'preprocessed', 'samples.fa')
    check_outfile(outfile)
    check_outfile(catted_samples)
    cat_samples(samples_file, catted_samples)

    # multiprocessing queue for transferring data to the main thread
    q = Manager().Queue()

    process_file(process, catted_samples, q=q, k=k)
    
    counter = Counter()
    while not q.empty():
        counter.update(q.get())
    printd('Finished consolidating counters.')
   
    # write counter to file
    write_dict(counter, outfile, sep='\t')

    # remove catted samples file
    remove(catted_samples)

if __name__ == '__main__':
    parse_args()
    main()

