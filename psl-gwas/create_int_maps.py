import pickle
from utility import get_params, parse_args, printd

def create_sample_int_map(samples, sim_file):
    printd('Creating sample int map...')
    sim = {}
    with open(samples, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines[1:]): # ignore header
        name = line.split('\t')[0]
        sim[name] = i
        sim[i] = name
    with open(sim_file, 'wb') as f:
        pickle.dump(sim, f)
    printd('Successfully created sample int map.')


def create_phenos_int_map(phenos, pim_file):
    printd('Creating pheno int map...')
    df = pd.read_csv(phenos, sep='\t')
    with open(phenos, 'r') as f:
        line = f.readline() # only need first line which contains headers
    pim = {}
    for i, p in enumerate(line[1:]): # skip id column
        pim[i] = p
        pim[p] = i
    with open(pim_file, 'wb') as f:
        pickle.dump(pim, f)
    printd('Successfully created phenos int map.')


def create_unitig_int_map(unitigs, uim_file):
    printd('Creating unitig int map...')
    uim = {}
    with open(unitigs, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        u = line.split()[0]
        uim[u] = i
        uim[i] = u
    with open(uim_file, 'wb') as f:
        pickle.dump(uim, f)
    printd('Successfully created unitig int map.')

def main():
    params = get_params()
    project = params['project']
    samples = params['samples']
    phenos = params['phenos']
    unitigs = f'{project}/data/preprocessed/unitig_sample_map.txt'
    sim_file = f'{project}/data/preprocessed/sample_int_map.pkl'
    pim_file = f'{project}/data/preprocessed/pheno_int_map.pkl'
    uim_file = f'{project}/data/preprocessed/unitig_int_map.pkl'

    create_sample_int_map(samples, sim_file)
    create_phenos_int_map(phenos, pim_file)
    create_unitig_int_map(unitigs, uim_file)

if __name__ == '__main__':
    parse_args()
    main()
