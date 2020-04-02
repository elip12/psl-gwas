#!/usr/bin/env python3
from multiprocessing import Pool
from os.path import getsize, isfile
import pickle
import argparse
import yaml

"""
Many helper methods for processing large files.
"""

DEBUG = False
PARAMS = None

    # parse args
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true',
        help='turn on debug mode')
    parser.add_argument('--project', required=True, type=str,
        help='name of project, defined with startproject.sh <name>')
    parser.add_argument('--samples', required=True, type=str,
        help='basename of samples file. ex: samples.tsv')
    parser.add_argument('--phenos', required=True, type=str,
        help='basename of phenos file. ex: phenos.tsv')
    parser.add_argument('--threads', default=2, type=int,
        help='max number of threads used concurrently')
    parser.add_argument('--mem', default=12, type=int,
        help='max amount of memory used concurrently (GB)')
    parser.add_argument('-k', '--k', default=30, type=int,
        help='kmer length in nucleotide bases')
    parser.add_argument('--upperfreq', default=0.98, type=float,
        help='kmer length in nucleotide base')
    parser.add_argument('--lowerfreq', default=0.02, type=float,
        help='kmer length in nucleotide bases')
    parser.add_argument('--thresh', default=5, type=int,
        help='kmer length in nucleotide bases')
    parser.add_argument('-p', '--param', action='store_true',
        help=('ignore threads, mem, k, lowerfreq, upperfreq, and thresh options'
            ' and use param file in project directory'))
    args = parser.parse_args()
    global DEBUG
    DEBUG = args.debug
    global PARAMS
    if args.param:
        PARAMS = read_yaml(args.param).update({'project': args.project})
    else:
        PARAMS = vars(args)

def read_yaml(fname):
    with open(fname, 'r') as f:
        return yaml.safe_load(f)

# prints only if debug mode is on
def printd(*args, **kwargs):
    if DEBUG == True:
        print(*args, **kwargs)

# loads a pickle file and returns the data
def load_pickle(fname):
    printd(f'Loading pickled data from {fname}')
    with open(fname, 'rb') as f:
        data = pickle.load(f)
    return data

# dumps data into a pickle file
def dump_pickle(data, fname):
    printd(f'Dumping pickled data to {fname}')
    with open(fname, 'wb') as f:
        pickle.dump(data, f)

# writes a data to a file, where each line in the file has one entry's
# data in the dictionary, formatted as {k}{v}. Optional separator
# to separate keys and values
def write_dict(data, fname, sep=''):
    printd(f'Writing data to {fname}')
    with open(fname, 'a+') as f:
        f.writelines(f'{k}{sep}{v}\n' for k,v in data.items())

# writes data to a file, where each line corresponds to an element in
# the data list
def write_list(data, fname):
    printd(f'Writing data to {fname}')
    with open(fname, 'a+') as f:
        s = '\n'.join(data) + '\n'
        f.write(s)

# reads a chunk of size chunk_size from fname into a string,
# splits it into a list, and calls process
def process_chunk(process_fn, fname, chunk_start, chunk_size, n,
        args, kwargs):
    printd(f'Starting chunk {n}')
    with open(fname, 'r') as f:
        f.seek(chunk_start)
        lines = f.read(chunk_size).splitlines()
    process_fn(lines, *args, **kwargs)

# breaks input file into chunks to minimize reads
def chunkify(fname, size=1024):
    file_end = getsize(fname)
    with open(fname,'rb') as f:
        chunk_end = f.tell()
        while chunk_end <= file_end:
            chunk_start = chunk_end
            f.seek(size, 1)
            f.readline()
            chunk_end = f.tell()
            chunk_size = chunk_end - chunk_start
            yield chunk_start, chunk_size

    # processes the file in chunks
def process_file(process_fn, fname, *args, **kwargs):
    input_size = getsize(fname)
    num_workers = PARAMS['threads']
    with Pool(num_workers) as pool:
        jobs = []
        n = 0

        optimal_size = int(input_size / num_workers) + 1
        max_mem_per_thread = int(10**9 * PARAMS['mem'] / PARAMS['threads'])
        size = min(optimal_size, max_mem_per_thread)
        
        for chunk_start, chunk_size in chunkify(fname, size):
            n += 1
            jobs.append(pool.apply_async(process_chunk,
                (process_fn, fname, chunk_start, chunk_size, n, args, kwargs)))

        # wait for all jobs to finish
        n = 0
        for job in jobs:
            job.get()
            n += 1
            printd(f'Finished chunk {n}')    

def get_params():
    return PARAMS

def file_exists(outfile):
    if isfile(outfile):
        printd(f'File {outfile} exists; skipping.')
        return 1
    return 0
