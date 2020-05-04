#!/usr/bin/env python3
###############################################################################
##  utility.py
##  This file holds many helper functions for processing large files.
##  Here is how it works:
##
##  Other scripts call parse_args(), defined here, to parse the command-line
##  arguments into a global PARAMS variable which they can then access through
##  the get_params() function, defined here. In this way, this file can access
##  command line args like the max number of threads allowed.
##  
##  Other scripts also define their own main() and process() functions.
##  In their main method, they invoke process_file(), defined here, passing
##  in their process() function, an input file, and some list and/or kwargs
##  corresponding to the parameters of their process function.
##
##  process_file() takes in, at minimum, a process function and an input file.
##  It splits the input file into chunks. The number of chunks is derived from
##  the number of allowed threads, and the max allowed memory. It then invokes
##  the script's process() function on that chunk. This allows other scripts
##  to efficiently utilize python's multiprocessing with very little code
##  or complexity.
###############################################################################
from multiprocessing import Pool
from os.path import getsize, isfile, join
import pickle
import argparse
import yaml
from random import randint

# global complement map, only needs to be created once
COMPLEMENT_MAP = str.maketrans('ATCG', 'TAGC')

DEBUG = False
PARAMS = None

# parse args and set global DEBUG and PARAMS vars
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--debug', action='store_true',
        help='turn on debug mode')
    parser.add_argument('--project', required=True, type=str,
        help='name of project, defined with startproject.sh <name>')
    parser.add_argument('--sample', required=True, type=str,
        help='basename of samples file. ex: samples.tsv')
    parser.add_argument('--pheno', required=True, type=str,
        help='basename of phenos file. ex: phenos.tsv')
    parser.add_argument('--threads', default=24, type=int,
        help='max number of threads used concurrently')
    parser.add_argument('--mem', default=350, type=int,
        help='upper memory limit (GB)')
    parser.add_argument('-k', '--k', default=30, type=int,
        help='kmer length in nucleotide bases')
    parser.add_argument('--minkf', default=0.01, type=float,
        help='minimum kmer frequency')
    parser.add_argument('-maxkf', default=0.95, type=float,
        help='maximum kmer frequency')
    parser.add_argument('--thresh', default=0.5, type=float,
        help='penetrance threshold for filtering in preprocessing')
    parser.add_argument('-p', '--param', action='store_true',
        help=('ignore sample, pheno, truth, threads, mem, k, minkf, maxkf,'
            ' and thresh options and use param file in project directory'))
    parser.add_argument('--truth', type=str,
        help=('fasta file holding truths data for benchmarking.'
            'Labels correspond to phenos, sequences hold genes or unitigs'
            'that cause the phenotype'))
    args = parser.parse_args()
    global DEBUG
    DEBUG = args.debug
    global PARAMS
    PARAMS = vars(args)
    if args.param:
        PARAMS.update(read_yaml(join(args.project, 'parameters.yaml')))

# read a yaml file and return the python dict
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

# writes a dict to a file, where each line in the file has one entry's
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

def write_2_files(data1, file1, data2, file2, lock):
    lock.acquire()
    myid = randint(1,100)
    printd('Lock acquired', myid)
    try:
        if file1 is not None:
           write_list(data1, file1)
        if file2 is not None:
            write_list(data2, file2)
    except Exception as e:
        print('Error: unable to write to files:', e)
    finally:
        lock.release()
        printd('Lock released', myid)

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
        max_mem_per_thread = int(10**9 / 2 * PARAMS['mem'] / PARAMS['threads'])
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

# allows other scripts to access PARAMS
def get_params():
    return PARAMS

# checks existence of a file and printds a warning if it does not
def file_exists(outfile):
    if isfile(outfile):
        printd(f'Found file {outfile}.')
        return 1
    return 0

# complements a kmer
def complement(kmer):
   return kmer.translate(COMPLEMENT_MAP)
