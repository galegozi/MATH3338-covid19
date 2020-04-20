# Problem: Find SARS-CoV-2.
import Levenshtein as L
import random as R
from multiprocessing import Pool
# Initial population: SARS+single mutation+padding

# alphabet = ['A', 'C', 'G', 'T']


def mutate(seq):
    if type(seq) != str:
        import sys
        print(seq)
        sys.exit("Trying to mutate on an object with a bad type")
    # a = alphabet
    alphabet = ['A', 'C', 'G', 'T']
    # choose a character to mutate
    pos = R.randrange(len(seq))
    if seq[pos] in alphabet:
        alphabet.remove(seq[pos])
    ch = R.choice(alphabet)
    # alphabet = ['A', 'C', 'G', 'T']
    return seq[:pos] + ch + seq[pos+1:]


def padding_worker_fxn(x):
    alphabet = ['A', 'C', 'G', 'T']
    return ''.join([R.choice(alphabet) for _ in range(x)])


def padding(length, seq, work_pool=None, workers=0):
    if type(seq) != str:
        print(seq)
        import sys
        sys.exit("Seq is not a string in padding")
    alphabet = ['A', 'C', 'G', 'T']
    if work_pool:
        return seq + ''.join(work_pool.map(padding_worker_fxn,
                                           [(length * (r+1))//workers - (length * r)//workers for r in range(workers)]))
    return seq + "".join([R.choice(alphabet) for _ in range(length)])

# Fitness is Levenshtine distance


def pop_fit_work_fxn(arg):
    (x, target) = arg
    if type(x) == str:
        print(x)
        import sys
        sys.exit("x is a string in pop_fit_work_fxn")
    if type(target) != str:
        print(target)
        import sys
        sys.exit("Target is not a string in pop_fit_work_fxn")
    return sum(L.distance(p, target) for p in x)


def pop_fitness(pop, target, work_pool=None, workers=0):
    for p in pop:
        if type(p) != str:
            print(p)
            import sys
            sys.exit("The population contains a non-string element in pop_fitness")
    if type(target) != str:
        print(target)
        import sys
        sys.exit("The target is not a string")
    if work_pool:
        l = len(pop)
        return sum(
            work_pool.map(
                pop_fit_work_fxn,
                [(pop[(l*r)//workers:(l*(r+1))//workers], target)
                 for r in range(workers)]
            )
        )/l
    return sum(L.distance(p, target) for p in pop)/len(pop)


def best_fit_helper(arg):
    (x, target) = arg
    if type(target) != str:
        print(target)
        import sys
        sys.exit("target is not a string in best_fit_helper")
    for k in x:
        if type(k) != str:
            print(k)
            import sys
            sys.exit("x contains a non-string argument in best_fit_helper")
    return min(L.distance(p, target) for p in x)


def best_fit(pop, target, work_pool=None, workers=0):
    for p in pop:
        if type(p) != str:
            import sys
            print(p)
            sys.exit("pop contains an item that is not a string in best_fit")
    if type(target) != str:
        import sys
        print(target)
        sys.exit("target is not a string in best_fit")
    if work_pool:
        l = len(pop)
        return min(
            work_pool.map(
                best_fit_helper,
                [(pop[(l*r)//workers:(l*(r+1))//workers], target)
                 for r in range(workers)]
            )
        )
    return min(L.distance(p, target) for p in pop)


def worst_helper(arg):
    (x, target) = arg
    for e in x:
        if type(e) != str:
            print(e)
            import sys
            sys.exit("x contains a bad element in worst_helper")
    if type(target) != str:
        print(target)
        import sys
        sys.exit("target is not a string in worst_helper")
    return max(L.distance(p, target) for p in x)


def worst_fit(pop, target, work_pool=None, workers=0):
    if work_pool:
        l = len(pop)
        return max(
            work_pool.map(
                worst_helper,
                [(pop[(l*r)//workers:(l*(r+1))//workers], target)
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


def distance_calc_worker(arg):
    (x, target) = arg
    return [(L.distance(p, target), p) for p in x]


def snd(x):
    return [e[1] for e in x]


def build_parent(my_list):
    output = []
    for ind in my_list:
        if 0.1 > R.random():
            output.append(ind)
    return output

def list_mutate(my_list):
    for i in range(len(my_list)):
        if type(my_list[i]) != str:
            print(my_list[i])
            import sys
            sys.exit("MyList[%d] is not a string in list_mutate" % i)
        if 0.1 > R.random():
            my_list[i] = mutate(my_list[i])
    return my_list

def gen_next_pop(pop, target, retain=0.2, random_select=0.1, mutate=0.1, work_pool=None, workers=0):
    # graded = [(L.distance(p, target), p) for p in pop]
    l = len(pop)
    if l != 1000:
        import sys
        sys.exit("Wrong population size in gen_next_pop")
    g = work_pool.map(
        distance_calc_worker,
        [(pop[(l*r)//workers:(l*(r+1))//workers], target)
         for r in range(workers)]
    )
    graded = [item for sub in g for item in sub]
    # for e in g:
    #     graded += e
    if len(graded) != 1000:
        import sys
        sys.exit("Population size changed during grading in gen_next_pop")
    graded.sort()
    # graded = [x[1] for x in sorted(graded)]
    g = work_pool.map(
        snd,
        [graded[(l*r)//workers:(l*(r+1))//workers] for r in range(workers)]
    )
    graded = []
    for e in g:
        graded += e
    if len(graded) != 1000:
        import sys
        sys.exit("Population size changed during snd call in gen_next_pop")
    keep = int(len(graded)*retain)
    parents = graded[:keep]
    a_list = graded[keep:]
    l = len(a_list)
    temp = work_pool.map(
        build_parent,
        [a_list[(l*r)//workers:(l*(r+1))//workers] for r in range(workers)]
    )
    for t in temp:
        parents += t
    # for ind in graded[keep:]:
    #     if random_select > R.random():
    #         parents.append(ind)
    # for i in range(len(parents)):
    #     if mutate > R.random():
    #         parents[i] = mutate(parents[i])
    l = len(parents)
    temp = work_pool.map(
        list_mutate,
        [parents[(l*r)//workers:(l*(r+1))//workers] for r in range(workers)]
    )
    parents = []
    for t in temp:
        parents += t
    ind_to_add = len(pop) - len(parents)
    children = []
    # TODO: Paralelize the following loop.
    temp = work_pool.map(
        children_builder,
        [(parents, (ind_to_add * (r+1))//workers - (ind_to_add * r)//workers) for r in range(workers)]
    )
    # while(len(children) < ind_to_add):
    #     p1 = R.choice(parents)
    #     p2 = R.choice(parents)
    #     if p1 != p2:
    #         child = breed(p1, p2, work_pool, workers)
    #         children.append(child)
    parents += children
    return parents

def children_builder(arg):
    (parents, count) = arg
    children = []
    for _ in range(count):
        p1 = R.choice(parents)
        p2 = R.choice(parents)
        while p1 == p2:
            p1 = R.choice(parents)
            p2 = R.choice(parents)
        children.append(breed(p1, p2))


if __name__ == "__main__":
    print('Done Imports')
    # First, read in files
    f = open("genomes/sars_cov_2/sars_cov_2_ref.fasta", "r")
    covid_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
    f.close()
    if type(covid_seq) != str:
        print(covid_seq)
        import sys
        sys.exit("Covid sequence")
    print('Done reading the COVID-CoV-2 sequence')
    f = open("genomes/sars_urbani/default_sars.fasta", "r")
    sars_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
    f.close()
    if type(sars_seq) != str:
        print(sars_seq)
        import sys
        sys.exit("Sars Sequence")
    print('Done reading the SARS Urbani sequence')
    # print(sars_seq)
    # print(covid_seq)
    # Start by applying 500 mutations on sars virus
    work_pool = Pool(4)
    gen = list(map(mutate, [sars_seq]*1000))  # TODO: Multi-thread this.
    for g in gen:
        if type(g) != str:
            print(gen)
            import sys
            sys.exit("Making gen")
    p = abs(len(covid_seq) - len(sars_seq))
    print("here is p: ", p)
    gen = [padding(p, seq) for seq in gen]
    for g in gen:
        if type(g) != str:
            print(gen)
            import sys
            sys.exit("adding padding")
    print("Before the simulation, gen has %d elements." % len(gen))
    # [g1, g2, g3, g4] = work_pool.map(gen_mutations, [(sars_seq, 125)]*4)
    # gen = g1+g2+g3+g4
#     print(gen)
    # Expect 6013 with current files
    min_dist = L.distance(covid_seq, sars_seq)  # If everything so far works, then this works too.
    f = open("output/simulation/sim.csv", "w")
    f.write("Closest, Current Closest\n")
    f.close()
    print('Ready for the simulation')
    for i in range(50000):
        gen = gen_next_pop(gen, covid_seq, work_pool=work_pool, workers=4)
        for g in gen:
            if type(g) != str:
                print(gen)
                import sys
                sys.exit("gen failed in generation %d (one of the elements is not a string)" % i)
        my_file = "output/backup/gen%d.txt" % i
        f = open(my_file, "w")
        for g in gen:
            f.write("%s\n" % g)
        f.close()
        if len(gen) != 1000:
            import sys
            sys.exit("Error: the length of gen is not right. Please check the log file for generation %d" % i)
        best = best_fit(gen, covid_seq, work_pool, 4)
        min_dist = min(best, min_dist)
        avg = pop_fitness(gen, covid_seq, work_pool, 4)
        worst = worst_fit(gen, covid_seq, work_pool, 4)
        f = open("output/simulation/sim.csv", "a")
        f.write("%d, %d, %f, %d" % (min_dist, best, avg, worst))
        f.close()
        print(i, min_dist, avg, best, worst)