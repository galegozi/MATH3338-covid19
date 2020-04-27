def lcs(X, Y):
    m = len(X)
    n = len(Y)
    L = [['' for i in range(n+1)] for j in range(2)]
    bi = bool

    for i in range(m):
        bi = i & 1

        for j in range(n+1):
            if (i == 0 or j == 0):
                L[bi][j] = ''

            elif (X[i] == Y[j - 1]):
                L[bi][j] = L[1 - bi][j - 1] + X[i]

            else:
                s1 = L[1 - bi][j]
                s2 = L[bi][j - 1]
                if len(s1) >= len(s2):
                    L[bi][j] = s1
                else:
                    L[bi][j] = s2
    return L[bi][n]
# load sequences from files.
f = open("genomes/sars_cov_2/sars_cov_2_ref_NC_045512.fasta", "r")
covid_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
f = open("genomes/RaTG13_MN996532.fasta", "r")
ratg13_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
f.close()
# Compute LCS on sequences.
seq = lcs(covid_seq, ratg13_seq)
# Fill in missing characters.
covid_pos = 0
rat_pos = 0
new_covid = ''
new_rat = ''
# for ch in seq:
#     print(ch)
f = open("lcs.txt", "w")
f.write(seq)
f.close()
print(len(seq))
