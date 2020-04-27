def lcs(X, Y):
    m = len(X)
    n = len(Y)
    L = [[None]*(n+1) for i in range(m+1)]
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i][j] = ''
            elif X[i-1] == Y[j-1]:
                L[i][j] = L[i-1][j-1]+X[i-1]
            else:
                s1 = L[i-1][j]
                s2 = L[i][j-1]
                if len(s1) >= len(s2):
                    L[i][j] = s1
                else:
                    L[i][j] = s2
    return L[m][n]

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
print(len(seq))