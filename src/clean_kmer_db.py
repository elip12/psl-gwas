import multiprocessing as mp, os, pickle, argparse

DEBUG = False

def printd(*args, **kwargs):
    if DEBUG == True:
        print(*args, **kwargs)

def read_cim():
    with open('../data/preprocess/class_int_map.pickle', 'rb') as f:
        cim = pickle.load(f)
    return cim

def write_kmer_class(lines):
    with open('data/psl/resistance_kmer_class.psl', 'a+') as f:
        s = '\n'/join(lines) + '\n'
        f.write(s)

def write_kmers(lines):
    printd('Writing kmers to file...')
    with open('data/output_cleaned.txt', 'a+') as f:
        s = '\n'.join(lines) + '\n'
        f.write(s)

def process(lines, cim):
    lines_to_write = []
    resistance_kmer_class_lines = []
    printd('Removing location metadata')
    for line in lines:
        linelist = line.split()
        if len(linelist) < 11:
            continue
        line_to_write = []
        for csv in linelist:
            values = csv.split(',')
            line_to_write.append(values[0])
        lines_to_write.append(' '.join(line_to_write))
        rkcl = [linelist[0] + ' ' cim[c] for c in cim]
        resistance_kmer_class_lines.append('\n'.join(rkcl))
    write_kmers(lines_to_write)
    printd('Done.')

# creates kmer dict from input chunk
def process_wrapper(chunkStart, chunkSize, cim): 
    with open('data/input.txt') as f:
        f.seek(chunkStart)
        lines = f.read(chunkSize).splitlines()
    process(lines, cim)

# breaks input file into chunks to minimize reads
def chunkify(fname, size=1024):
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

    NUM_WORKERS = 16

    #initialize raw data, multiprocessing pool, and jobs queue
    pool = mp.Pool(NUM_WORKERS)
    jobs = []

    cim = read_cim

    # create jobs
    n = 0
    for chunkStart,chunkSize in chunkify('data/input.txt', int(0.8 * 10**9)):
        n += 1
        printd(f'Starting chunk {n}...')
        jobs.append(pool.apply_async(process_wrapper, (chunkStart,chunkSize, cim)))

    # wait for all jobs to finish
    n = 0
    for job in jobs:
        job.get()
        n += 1
        printd(f'Finished chunk {n}...')

    pool.close()
