# get ref genome from genomes folder.
# for each genome, compute the levenstine distance to the reference.
# output csv with id, lev(id)
import Levenshtein as L
import pandas as pd
f = open("genomes/sars_cov_2_ref.fasta", "r")
ref_seq = ''.join([x[:-1] for x in f.readlines()[1:]])
# dist = Levenshtein.distance(ref_seq, 'abc')
f.close()
f = open("genomes/sars_cov_2.fasta")
genomes = f.read().split('>')[1:]
f.close()
g = pd.read_csv("genomes/sars_cov_2.csv")
merge = g[['Accession', 'Collection_Date', 'Nuc_Completeness']]
merge.drop(merge[merge['Nuc_Completeness'] != 'complete'].index, inplace=True)
merge['Distance'] = [0]*len(merge.index)
print(merge)
# dataset = {
#     'Accession': [],
#     'Distance': []
# }
# print(merge[merge['Accession'] == 'MT308702'])
# print("loading", end='')
count = 0
for g in genomes:
    code = g.split()[0]
    if(code not in list(merge['Accession'])):
        continue
    lines = g.split('\n')
    info = lines[0]
    genetic = ''.join(lines[1:])
    distance = L.distance(ref_seq, genetic)
    index = merge[merge['Accession'] == code].index.tolist()[0]
    # merge[merge['Accession'] == code]['Distance'] = distance
    merge.at[index, 'Distance'] = distance
    print('count = ', count)
    print(merge)
    count += 1
    # print(code, L.distance(ref_seq, genetic))
    # merge['Distance'] = distance
# print(dataset)
print(merge)
merge.to_csv('output/sars_cov_2.csv', index = False)
# print(dist)