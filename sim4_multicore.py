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
    # Choose a character to replace seq[pos] with.
    ch = R.choice(alphabet)
    # Return the string with the substitution.
    return seq[:pos] + ch + seq[pos+1:]


def padding(x):
    # Creating a random string of length x with characters from the alphabet.
    alphabet = ['A', 'C', 'G', 'T']
    return ''.join([R.choice(alphabet) for _ in range(x)])

def run_me(arg):
    # Because multiprocessing works with one argument.
    (covid_seq, ratg13_seq) = arg
    # How many characters of padding do we need?
    num_padding = len(covid_seq) - len(ratg13_seq)
    # Add the padding.
    ratg13_seq += padding(num_padding)
    # Whats the current distance?
    dist_now = L.distance(covid_seq, ratg13_seq)
    # Here is the current sequence we are working with.
    current_seq = ratg13_seq
    # We start with the RATG13 sequence.
    worsening = 0
    # Counting the nuber of iterations where the mutations tried to increase the distance.
    for iter in range(50000):
        # Go through 50000 iterations
        if iter % 100 == 0:
            # Printing iteration counts for analysis/debug purposes.
            print("Iteration %d" % iter)
        mutated = mutate(current_seq)
        # Apply a mutation to the sequence.
        mutated_dist = L.distance(mutated, covid_seq)
        # Compute the distance after the mutation was applied.
        if mutated_dist <= dist_now:
            # Has our distance decreased?
            # If so, fix the distance and the sequence.
            dist_now = mutated_dist
            current_seq = mutated
        else:
            # We have gotten worse.
            worsening += 1
        print("The distance is now %d" % dist_now)
    print("The final distance is %d units away from COVID 19." % dist_now)
    print("In %d iterations, the distance tried to get worse." % worsening)

if __name__ == "__main__":
    # Open the files for the sequences and the LCS
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
    p = Pool(4)
    p.map(run_me, [(covid_seq, new_rat)]*4)