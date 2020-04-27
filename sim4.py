import random as R


def padding(x):
    alphabet = ['A', 'C', 'G', 'T']
    return ''.join([R.choice(alphabet) for _ in range(x)])


# Load sequences and LCS from file
f = open("genomes/sars_cov_2/sars_cov_2_ref_NC_045512.fasta", "r")
covid_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
f = open("genomes/RaTG13_MN996532.fasta", "r")
ratg13_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
f = open('lcs.txt', 'r')
lcs = f.read()
# Fill in missing characters.
covid_pos = 0
rat_pos = 0
new_rat = ''
for ch in lcs:
    while covid_seq[covid_pos] != ch and ratg13_seq[rat_pos] != ch:
        new_rat += ratg13_seq[rat_pos]
        covid_pos += 1
        rat_pos += 1
    if covid_seq[covid_pos] == ch and ratg13_seq[rat_pos] == ch:
        new_rat += ch
        covid_pos += 1
        rat_pos += 1
        continue
    while covid_seq[covid_pos] != ch:
        new_rat += covid_seq[covid_pos]
        covid_pos += 1
    while ratg13_seq[rat_pos] != ch:
        # new_rat += R.choice(['A', 'C', 'G', 'T'])
        rat_pos += 1
    new_rat += ch
    covid_pos += 1
    rat_pos += 1
# Pad new_rat
new_rat += padding(len(covid_seq) - len(new_rat))
#then mutate.