from random import randrange, choice
import Levenshtein as L
# First, read in files.
f = open("genomes/sars_cov_2/sars_cov_2_ref.fasta", "r")
covid_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
f = open("genomes/sars_urbani/default_sars.fasta", "r")
sars_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
# print(sars_seq)
# print(covid_seq)
# Start by applying 500 mutations on sars virus
def mutate(seq):
    # choose a character to mutate
    pos = randrange(len(seq))
    # What do I do with that position?
    todo = choice([
        "Remove",
        "Replace",
        "Add"
    ])
    if todo == "Remove":
        return seq[:pos] + seq[pos+1:]
    # either a replace or an add. Would need a character.
    l = ['A', 'C', 'G', 'T']
    l.remove(seq[pos])
    ch = choice(l)
    if(todo == "Replace"):
        return seq[:pos] + ch + seq[pos+1:]
    return seq[:pos] + ch + seq[pos:]
gen = []
for _ in range(500):
    gen.append(mutate(sars_seq))
for _ in range(50000):
    gen = list(map(lambda x: (L.distance(sars_seq, x), x), gen))
    gen.sort()
    gen.map(lambda x: x[1], gen)
    remutate = []
    
    print(gen[0])
# for _ in range(1000):
#     print(L.distance(sars_seq, mutate(sars_seq)))