from random import randrange, choice
import Levenshtein as L
from multiprocessing import Pool

def gen_mutations(arg):
    (seq, l) = arg
    return list(map(lambda x: mutate(seq), [1]*l))

# def mutate(seq):
#     # choose a character to mutate
#     pos = randrange(len(seq))
#     # What do I do with that position?
#     todo = choice([
#         "Remove",
#         "Replace",
#         "Add"
#     ])
#     if todo == "Remove":
#         return seq[:pos] + seq[pos+1:]
#     # either a replace or an add. Would need a character.
#     l = ['A', 'C', 'G', 'T']
#     if seq[pos] in l:
#         l.remove(seq[pos])
#     ch = choice(l)
#     if(todo == "Replace"):
#         return seq[:pos] + ch + seq[pos+1:]
#     return seq[:pos] + ch + seq[pos:]

def mutate(seq):
    # choose a character to mutate
    pos = randrange(len(seq))
    # What do I do with that position?
    l = ['A', 'C', 'G', 'T']
    if seq[pos] in l:
        l.remove(seq[pos])
    ch = choice(l)
    return seq[:pos] + ch + seq[pos+1:]


def get_dist(arg):
    (l, covid_seq) = arg
    return list(map(lambda x: (L.distance(covid_seq, x), x), l))


def getFElement(l):
    return list(map(lambda x: x[0], l))


def getSElement(l):
    return list(map(lambda x: x[1], l))


def list_mutate(l):
    return list(map(lambda k: mutate(k), l))


if __name__ == "__main__":
    print('Done Imports')
    # First, read in files.
    f = open("genomes/sars_cov_2/sars_cov_2_ref.fasta", "r")
    covid_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
    f.close()
    print('Done reading the COVID-CoV-2 sequence')
    f = open("genomes/sars_urbani/default_sars.fasta", "r")
    sars_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
    f.close()
    print('Done reading the SARS Urbani sequence')
    # print(sars_seq)
    # print(covid_seq)
    # Start by applying 500 mutations on sars virus
    work_pool = Pool(4)
    [g1, g2, g3, g4] = work_pool.map(gen_mutations, [(sars_seq, 125)]*4)
    gen = g1+g2+g3+g4
    min_dist = L.distance(covid_seq, sars_seq) # Expect 6013 with current files
    f = open("output/simulation/sim.csv", "w")
    f.write("Closest, Current Closest\n")
    f.close()
    print('Ready for the simulation')
    for _ in range(50000):
        gen_len = len(gen)
        g1 = gen[:gen_len//4]
        g2 = gen[gen_len//4:gen_len//2]
        g3 = gen[gen_len//2:3*gen_len//4]
        g4 = gen[3*gen_len//4:]
        [g1, g2, g3, g4] = work_pool.map(get_dist, [(g1, covid_seq), (g2, covid_seq), (g3, covid_seq), (g4, covid_seq)])
        # gen = g1 + g2 + g3 + g4
        # gen = list(map(lambda x: (L.distance(covid_seq, x), x), gen))
        [d1, d2, d3, d4] = work_pool.map(getFElement, [g1, g2, g3, g4])
        # distances = list(map(lambda x: x[0], gen))
        [m1, m2, m3, m4] = work_pool.map(min, [d1, d2, d3, d4])
        current_min = min(m1, m2, m3, m4)
        min_dist = min(min_dist, m1, m2, m3, m4)
        gen = g1+g2+g3+g4
        gen.sort()
        g1 = gen[:gen_len//4]
        g2 = gen[gen_len//4:gen_len//2]
        g3 = gen[gen_len//2:3*gen_len//4]
        g4 = gen[3*gen_len//4:]
        [g1, g2, g3, g4] = work_pool.map(getSElement, [g1, g2, g3, g4])
        # gen = list(map(lambda x: x[1], gen))
        remutate = []
        direct_copy = []
        gen = g1+g2+g3+g4
        for _ in range(250):
            chosen = choice(gen)
            gen.remove(chosen)
            remutate.append(chosen)
            direct_copy.append(gen[0])
            gen = gen[1:]
        mut_len = len(remutate)
        [m1, m2, m3, m4] = [remutate[:mut_len//4], remutate[mut_len//4:mut_len//2],
                            remutate[mut_len//2:3*mut_len//4], remutate[3*mut_len//4:]]
        [m1, m2, m3, m4] = work_pool.map(list_mutate, [m1, m2, m3, m4])
        # remutate = list(map(lambda x: mutate(x), remutate))
        gen = remutate + direct_copy
        g1 = gen[:gen_len//4]
        g2 = gen[gen_len//4:gen_len//2]
        g3 = gen[gen_len//2:3*gen_len//4]
        g4 = gen[3*gen_len//4:]
        [g1, g2, g3, g4] = work_pool.map(list_mutate, [g1, g2, g3, g4])
        gen = g1 + g2 + g3 + g4
        f = open("output/simulation/sim.csv", "a")
        f.write("%d,%d\n" % (min_dist, current_min))
        f.close()
        print(min_dist, current_min)
    # for _ in range(1000):
    #     print(L.distance(sars_seq, mutate(sars_seq)))