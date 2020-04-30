import random as R
M=30000
S=30000
A=4
mutations = []
for _ in range(10**6):
    mutations.append(M)
    if R.random() <= M/(S*(A-1)):
        M -= 1
f = open("mysim.txt", 'w')
for e in mutations:
    f.write(str(e) + "\n")
f.close()