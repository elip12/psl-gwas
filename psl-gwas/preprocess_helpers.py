
COMPLEMENT_MAP = str.maketrans('ATCG', 'TAGC')


def kmer_frequency(val, upper, lower):
    if lower <= x <= upper:
        return False
    return True

# for small files
def num_samples(fname):
    with open(fname, 'r') as f:                                                
        return len(f.readlines()) - 1

def parse_input(infile):
    seqs = {}
    with open(infile, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        sample, fname = line.rstrip().split('\t')
        if not fname.endswith('.fa'):
            continue
        with open(fname, 'r') as f:
            contiglines = f.readlines()
        seqs[sample] = [contigline.rstrip() for contigline in contiglines \
            if contigline[0] in ['A','T','C','G']]
   return seqs 

# DNA base complement
def complement(kmer):
   return kmer.translate(COMPLEMENT_MAP)

def consolidate(data, k):
    prev_linelist = data[0]
    prev_unitig = prev_linelist[0]
    prev_samples = prev_linelist[1]
    unitig_db_chunk = []
    for linelist in data[1:]:
        this_unitig = linelist[0]
        this_samples = linelist[1]
        
        # if kmers are sequential and the same set of samples contain both kmers
        if prev_unitig[-(k - 1):] == line[0:k - 1] \
                and len(this_samples) == len(prev_samples) \
                and set(this_samples) == set(prev_samples):
            this_unitig = f'{prev_unitig}{line[k - 1]}'
            prev_line = (this_unitig, prev_line[1])
        else:
            unitig_db_chunk.append(prev_line)
            prev_line = linelist
        
        prev_unitig = this_unitig
        prev_samples = this_samples
    return unitig_db_chunk

def filter_unitigs(data, THRESH, dfdisp, dfnodisp):
    kmer_db_chunk = []
    nphenos = dfdisp.shape[1]
    for linelist in data:
        kmer_db_line = [linelist[0]]
        disp = np.zeros(nphenos)
        nodisp = np.zeros(nphenos)

        for sample_id, _ in linelist[1]:
            if sample_id not in dfdisp.index:
                continue
            kmer_db_line.append(sample_id)
            # collect resistant/vulnerable frequencies for each antibiotic for
            # this kmer
            disp += dfdisp.loc[sample_id].to_numpy()
            nodisp += dfnodisp.loc[sample_id].to_numpy()

        # 1 test per antibiotic; kmer needs to pass only 1 to avoid
        # getting filtered out
        a = np.where((disp + nodisp >= THRESH) \
                    & (disp > 0) \
                    & (nodisp / disp < 0.05))[0]

        if a.size > 0:
            kmer_db_chunk.append('\t'.join(kmer_db_line))
    return kmer_db_chunk


# Creates a dictionary of all kmers passed in data, and their complements
# Then, iterates through genomes. If a genome contains a kmer, that genome's
# metadata is added to dict entry for that kmer.
# After all samples have been iterated over, writes kmers to file.
def create_unitig_sample_map(data, raw, k, q, upper, lower, thresh, dfdisp, dfnodisp):
    # get all kmers in chunk and complement them
    kmers = {}
    for line in data:
        kmer, count = line.split('\t')
        if kmer_frequency_fails(count, upper, lower):
            continue
        comp = complement_kmer(kmer)
        kmers[kmer] = []
        kmers[comp] = []
    
    # map all kmers in chunk to samples containing them
    for count, (raw_id, seq) in enumerate(raw.items()):
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= k: # ensure this contig is long enough to sample
                for i in range(l - k + 1):
                    kmer = contig[i: i + k]
                    if kmer in kmers:
                        kmers[kmer].append((raw_id, c_id))
    unitigs = consolidate(kmers.items(), k)
    unitigs = filter_unitigs(unitigs, thresh, dfdisp, dfnodisp)
    q.put(unitigs)
