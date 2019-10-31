from multiprocessing import Pool
from os.path import getsize
import pickle
import argparse

''' Description
A class that contains many helper methods for processing large files.
Subclasses need to implement the `process` method, and override main.
'''

class LargeFileProcessor():
    def __init__(self, num_workers, fname):
        self.num_workers = num_workers
        self.fname = fname
        try:
            self.input_size = getsize(fname)
        except Exception as e:
            print(f'{fname} is not a valid file')
            self.input_size = None
        
        # parse args
        self.parser = argparse.ArgumentParser()
        parser.add_argument('-d', action='store_true',
            help='turn on debug mode')
        args = self.parser.parse_args()
        self.DEBUG = args.d

    # prints only if debug mode is on
    def printd(self, *args, **kwargs):
        if self.DEBUG == True:
            print(*args, **kwargs)

    # loads a pickle file and returns the data
    def load_pickle(self, fname):
        self.printd(f'Loading pickled data from {fname}')
        with open(fname, 'rb') as f:
            data = pickle.load(f)
        return data

    # dumps data into a pickle file
    def dump_pickle(self, data, fname):
        self.printd(f'Dumping pickled data to {fname}')
        with open(fname, 'wb') as f:
            pickle.dump(data, f)

    # writes a data to a file, where each line in the file has one entry's
    # data in the dictionary, formatted as {k}{v}
    def write_dict(self, data, fname):
        self.printd(f'Writing data to {fname}')
        with open(fname, 'a+') as f:
            s = '\n'.join([f'{k}{v}' for k,v in data.items()]) + '\n'
            f.write(s)

    # writes data to a file, where each line corresponds to an element in
    # the data list
    def write_list(self, data, fname):
        self.printd(f'Writing data to {fname}')
        with open(fname, 'a+') as f:
            s = '\n'.join(data) + '\n'
            f.write(s)

    # Processes a single chunk of data, passed as a list of lines
    def process(self, data, *args, **kwargs):
        raise NotImplementedError

    # reads a chunk of size chunk_size from fname into a string,
    # splits it into a list, and calls process
    def process_chunk(self, fname, chunk_start, chunk_size, n, *args, **kwargs): 
        self.printd(f'Starting chunk {n}')
        with open(fname, 'r') as f:
            f.seek(chunk_start)
            lines = f.read(chunk_size).splitlines()
        self.process(lines, *args, **kwargs)

    # breaks input file into chunks to minimize reads
    def chunkify(self, fname, size=1024):
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
    def main(self, *args, **kwargs):
        jobs = []
        with Pool(self.num_workers) as pool:
            n = 0
            size = int(self.input_size / self.num_workers) + 1
            for chunk_start, chunk_size in self.chunkify(self.fname, size):
                n += 1
                jobs.append(pool.apply_async(process_chunk,
                    (chunk_start, chunk_size, n, *args, **kwargs)))

            # wait for all jobs to finish
            n = 0
            for job in jobs:
                job.get()
                n += 1
                printd(f'Finished chunk {n}')
    
    def main_wrapper(self):
        raise NotImplementedError

