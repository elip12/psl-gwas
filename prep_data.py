import multiprocessing as mp, os, pickle, argparse

DEBUG = False
pickle_file = 'data/raw.pickle'
K = 30

def printd(*args, **kwargs):
    if DEBUG == True:
        print(*args, **kwargs)

def load_raw():
    printd('Loading data...')
    with open(pickle_file, 'rb') as f:
        raw = pickle.load(f)
    return raw

def complement(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'G':
        return 'C'
    else:
        return 'G'

def write_kmers(kmers):
    printd('Writing kmers to file...')
    with open('data/output.txt', 'a+') as f:
        s = '\n'.join([f'{k}{v}' for k,v in kmers.items()]) + '\n'
        f.write(s)

def process(kmers):
    printd('Extracting kmer locations...')
    for count, (raw_id, seq) in enumerate(raw.items()):
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= K: # ensure this contig is long enough to sample
                for i in range(l - K + 1):
                    kmer = contig[i:i + K]
                    if kmer in kmers:
                        kmers[kmer] += f' {raw_id},{c_id},{i}'
        printd(f'\tProcessed genome {count + 1}', end='\r')
    write_kmers(kmers)
    printd('Done.')

# creates kmer dict from input chunk
def process_wrapper(chunkStart, chunkSize): 
    kmers = {}
    with open('data/input.txt') as f:
        f.seek(chunkStart)
        lines = f.read(chunkSize).splitlines()
        for line in lines:
            kmer = line.split(' ')[0]
            comp = ''.join(map(complement, kmer.split()))
            kmers[kmer] = ''
            kmers[comp] = ''
    process(kmers)
    # print(len(kmers))

# breaks input file into chunks to minimize reads
def chunkify(fname, size=150 * 10**6):
    fileEnd = os.path.getsize(fname)
    with open(fname,'rb') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs GWAS on microbial genomes')
    parser.add_argument('-d', action='store_true',
        help='turn on debug mode')
    args = parser.parse_args()
    DEBUG = args.d

    #initialize raw data, multiprocessing pool, and jobs queue
    raw = load_raw()
    #print(list(raw.items())[0])
    pool = mp.Pool(4)
    jobs = []

    # create jobs
    n = 0
    for chunkStart,chunkSize in chunkify('data/input.txt'):
        n += 1
       # if n > 4:
        #    break
        printd(f'Starting chunk {n}...')
        jobs.append(pool.apply_async(process_wrapper,(chunkStart,chunkSize)))

    # wait for all jobs to finish
    n = 0
    for job in jobs:
        job.get()
        n += 1
        printd(f'Finished chunk {n}...')

    pool.close()
