import pickle
def create_sample_int_map():
    sim = {}
    with open('data/sample_int_map.txt', 'r') as f:
        for line in f:
            i, s = line.split()
            s = s.split('.')[0]
            sim[s] = i
            sim[i] = s
    with open('data/sample_int_map.pickle', 'wb+') as f:
        pickle.dump(sim, f)
create_sample_int_map()
