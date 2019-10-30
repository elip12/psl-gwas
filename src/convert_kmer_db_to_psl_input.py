import multiprocessing as mp, os, pickle, argparse

DEBUG = False

def read_kmer_map():
    pass

# sample int map
def read_sim():
    with open('data/sample_int_map.pickle', 'rb') as f:
        sim = pickle.load(f)
    return sim

def printd(*args, **kwargs):
    if DEBUG == True:
        print(*args, **kwargs)

def write_lines(lines):
    with open('data/contains_sample_kmer.txt', 'a+') as f:
        s = '\n'.join(lines) + '\n'
        f.write(s)

def process(lines, sim):
    printd('Starting chunk')
    lines_to_write = []
    for line in lines:
        linelist = line.split()
        kmer = linelist[0]
        line_to_write = []
        for sample in linelist[1:]:
            line_to_write.append(kmer + ' ' + sim[sample] + ' 1.0')
        lines_to_write.append('\n'.join(line_to_write))
    write_lines(lines_to_write)

# creates kmer dict from input chunk
def process_wrapper(chunkStart, chunkSize, sim): 
    with open('data/input.txt') as f:
        f.seek(chunkStart)
        lines = f.read(chunkSize).splitlines()
    process(lines, sim)

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
    size = 1024
    #size = int(10**9 * 6.9 / NUM_WORKERS) + 1

    #initialize raw data, multiprocessing pool, and jobs queue
    pool = mp.Pool(NUM_WORKERS)
    jobs = []

    # read sample int map
    sim = read_sim()

    # create jobs
    n = 0
    for chunkStart,chunkSize in chunkify('data/input.txt', size):
        n += 1
        jobs.append(pool.apply_async(process_wrapper, (chunkStart, chunkSize, sim)))

    # wait for all jobs to finish
    n = 0
    for job in jobs:
        job.get()
        n += 1
        printd(f'Finished chunk {n}...')

    pool.close()
