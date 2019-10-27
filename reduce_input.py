with open('data/input_orig.txt', 'r') as fin:
    with open('data/input_reduced.txt', 'a+') as fout:
        for line in fin:
            s, n = line.split()
            if int(n) >= 10:
                fout.write(line)

