from multiprocessing import Pool
from os.path import getsize
import pickle
import argparse

"""
Many helper methods for processing large files.
"""

DEBUG = False

    # parse args
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', action='store_true',
        help='turn on debug mode')
    args = parser.parse_args()
    global DEBUG
    DEBUG = args.d

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
# data in the dictionary, formatted as {k}{v}
def write_dict(data, fname):
    printd(f'Writing data to {fname}')
    with open(fname, 'a+') as f:
        s = '\n'.join([f'{k}{v}' for k,v in data.items()]) + '\n'
        f.write(s)

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
def process_file(process_fn, num_workers, fname, *args, **kwargs):
    input_size = getsize(fname)
    with Pool(num_workers) as pool:
        jobs = []
        n = 0
        size = int(input_size / num_workers) + 1
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

