import pickle

def create_sample_int_map():
    sim = {}
    with open('data/intermediate/sample_int_map.txt', 'r') as f:
        for line in f:
            i, s = line.split()
            s = s.split('.')[0]
            sim[s] = i
            sim[i] = s
    with open('data/intermediate/sample_int_map.pickle', 'wb+') as f:
        pickle.dump(sim, f)

if __name__ == '__main__':
    create_sample_int_map()

