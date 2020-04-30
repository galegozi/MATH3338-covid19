import Levenshtein as L
import random as R
from multiprocessing import Pool


def padding(x):
    # Creating a random string of length x with characters from the alphabet.
    alphabet = ['A', 'C', 'G', 'T']
    return ''.join([R.choice(alphabet) for _ in range(x)])


f = open("genomes/sars_cov_2/sars_cov_2_ref_NC_045512.fasta", "r")
covid_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
f = open("genomes/RaTG13_MN996532.fasta", "r")
ratg13_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
f = open('lcs.txt', 'r')
lcs = f.read()
# Alignment
covid_pos = 0
rat_pos = 0
new_rat = ''
for ch in lcs:
    while covid_seq[covid_pos] != ch and ratg13_seq[rat_pos] != ch:
        # Do we have a sequence where neither COVID or RATG matches ch?
        # We retain the RATG version in the RATG sequence.
        new_rat += ratg13_seq[rat_pos]
        covid_pos += 1
        rat_pos += 1
    if covid_seq[covid_pos] == ch and ratg13_seq[rat_pos] == ch:
        # Have we reached a point where both the current character in COVID and RATG are currently ch from the LCS?
        new_rat += ch
        covid_pos += 1
        rat_pos += 1
        continue  # We skip the rest of the code for the current iteration.
    while covid_seq[covid_pos] != ch:
        # Does the COVID sequence exclusively have characters deviating from the LCS?
        # Then we fix it by adding the COVID characters to the RATG sequence (assume an insert)
        new_rat += covid_seq[covid_pos]
        covid_pos += 1
    while ratg13_seq[rat_pos] != ch:
        # We delete inserts exclusive to RATG.
        rat_pos += 1
    new_rat += ch
    covid_pos += 1
    rat_pos += 1
num_padding = len(covid_seq) - len(new_rat)
my_start = 0
# Add the padding.
iter = 10000
for _ in range(iter):
    # Whats the current distance?
    dist_now = L.distance(covid_seq, new_rat+padding(num_padding))
    print(dist_now)
    my_start += dist_now/iter
print("Average Starting Distance: %d" % my_start)
