import Levenshtein as L
import random as R
from multiprocessing import Pool


def mutate(seq):
    # a = alphabet
    alphabet = ['A', 'C', 'G', 'T']
    # choose a character to mutate
    pos = R.randrange(len(seq))
    if seq[pos] in alphabet:
        alphabet.remove(seq[pos])
    ch = R.choice(alphabet)
    return seq[:pos] + ch + seq[pos+1:]


def padding(x):
    alphabet = ['A', 'C', 'G', 'T']
    return ''.join([R.choice(alphabet) for _ in range(x)])

def run_me(arg):
    (covid_seq, ratg13_seq) = arg
    num_padding = len(covid_seq) - len(ratg13_seq)
    dist_now = L.distance(covid_seq, ratg13_seq)
    ratg13_seq += padding(num_padding)
    current_seq = ratg13_seq
    worsening = 0
    for _ in range(50000):
        mutated = mutate(current_seq)
        mutated_dist = L.distance(mutated, covid_seq)
        if mutated_dist <= dist_now:
            dist_now = mutated_dist
            current_seq = mutated
        if mutated_dist > dist_now:
            worsening += 1
        print("The distance is now %d" % dist_now)
    print("The final distance is %d units away from COVID 19." % dist_now)
    print("In %d iterations, the distance tried to get worse." % worsening)

if __name__ == "__main__":
    f = open("genomes/sars_cov_2/sars_cov_2_ref_NC_045512.fasta", "r")
    covid_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
    f.close()
    f = open("genomes/RaTG13_MN996532.fasta", "r")
    ratg13_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
    f.close()
    p = Pool(4)
    p.map(run_me, [(covid_seq, ratg13_seq)]*4)