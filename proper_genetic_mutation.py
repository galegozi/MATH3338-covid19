# Problem: Find SARS-CoV-2.
import Levenshtein as L
import random as R
from multiprocessing.pool import ThreadPool
# # Initial population: SARS+single mutation+padding


def mutate(seq):
    alphabet = ['A', 'C', 'G', 'T']
    # choose a character to mutate
    pos = R.randrange(len(seq))
    if seq[pos] in alphabet:
        alphabet.remove(seq[pos])
    ch = R.choice(alphabet)
    return seq[:pos] + ch + seq[pos+1:]


def padding(length, seq, work_pool=None, workers=0):
    alphabet = ['A', 'C', 'G', 'T']
    if work_pool:
        return seq + ''.join(work_pool.map(lambda x: ''.join([R.choice(alphabet) for _ in range(x)]),
                                           [(length * (r+1))//workers - (length * r)//workers for r in range(workers)]))
    return seq + "".join([R.choice(alphabet) for _ in range(length)])

# Fitness is Levenshtine distance


def pop_fitness(pop, target, work_pool=None, workers=0):
    if work_pool:
        l = len(pop)
        return sum(
            work_pool.map(
                lambda x: sum(L.distance(p, target) for p in x),
                [pop[(l*r)//workers:(l*(r+1))//workers]
                 for r in range(workers)]
            )
        )/l
    return sum(L.distance(p, target) for p in pop)/len(pop)


def best_fit(pop, target, work_pool=None, workers=0):
    if work_pool:
        l = len(pop)
        return min(
            work_pool.map(
                lambda x: min(L.distance(p, target) for p in x),
                [pop[(l*r)//workers:(l*(r+1))//workers]
                 for r in range(workers)]
            )
        )
    return min(L.distance(p, target) for p in pop)


def worst_fit(pop, target, work_pool=None, workers=0):
    if work_pool:
        l = len(pop)
        return max(
            work_pool.map(
                lambda x: max(L.distance(p, target) for p in x),
                [pop[(l*r)//workers:(l*(r+1))//workers]
                 for r in range(workers)]
            )
        )
    return max(L.distance(p, target) for p in pop)
# Population size 1000


def breed(p1, p2, work_pool=None, workers=0):
    if work_pool:
        l = len(p1)
        return ''.join(
            work_pool.map(single_breed, [
                (p1[(l*r)//workers:(l*(r+1))//workers], p2[(l*r)//workers:(l*(r+1))//workers]) for r in range(workers)
            ])
        )
    return ''.join(R.choice([p1[pos], p2[pos]]) for pos in range(len(p1)))


def single_breed(parents):
    (p1, p2) = parents
    return ''.join(R.choice([p1[pos], p2[pos]]) for pos in range(len(p1)))


def gen_next_pop(pop, target, retain=0.2, random_select=0.1, mutate=0.1, work_pool=None, workers=0):
    # graded = [(L.distance(p, target), p) for p in pop]
    l = len(pop)
    g = work_pool.map(
        lambda x: [(L.distance(p, target), p) for p in x],
        [pop[(l*r)//workers:(l*(r+1))//workers] for r in range(workers)]
    )
    graded = []
    for e in g:
        graded += e
    graded.sort()
    # graded = [x[1] for x in sorted(graded)]
    g = work_pool.map(
        lambda x: [e[1] for e in x],
        [graded[(l*r)//workers:(l*(r+1))//workers] for r in range(workers)]
    )
    keep = int(len(graded)*retain)
    parents = graded[:keep]
    # TODO: paralell processing for two following loops.
    for ind in graded[keep:]:
        if random_select > R.random():
            parents.append(ind)
    for i in range(len(parents)):
        if mutate > R.random():
            parents[i] = mutate(parents[i])
    ind_to_add = len(pop) - len(parents)
    children = []
    # TODO: Paralelize the following loop.
    while(len(children) < ind_to_add):
        p1 = R.choice(parents)
        p2 = R.choice(parents)
        if p1 != p2:
            child = breed(p1, p2, work_pool, workers)
            children.append(child)
    parents += children
    return parents


if __name__ == "__main__":
    print('Done Imports')
    # First, read in files
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
    work_pool = ThreadPool(4)
    gen = list(map(mutate, [sars_seq]*1000))
    p = abs(len(covid_seq) - len(sars_seq))
    gen = [padding(p, seq) for seq in gen]
    # [g1, g2, g3, g4] = work_pool.map(gen_mutations, [(sars_seq, 125)]*4)
    # gen = g1+g2+g3+g4
#     print(gen)
    # Expect 6013 with current files
    min_dist = L.distance(covid_seq, sars_seq)
    f = open("output/simulation/sim.csv", "w")
    f.write("Closest, Current Closest\n")
    f.close()
    print('Ready for the simulation')
    for _ in range(50000):
        gen = gen_next_pop(gen, covid_seq, work_pool=work_pool, workers=4)
        best = best_fit(gen, covid_seq, work_pool, 4)
        min_dist = min(best, min_dist)
        avg = pop_fitness(gen, covid_seq, work_pool, 4)
        worst = worst_fit(gen, covid_seq, work_pool, 4)
        f = open("output/simulations/sim.csv", "a")
        f.write("%d, %d, %f, %d" % (min_dist, best, avg, worst))
        print(avg,best, worst)
