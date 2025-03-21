import random
import sys

genelist=open(sys.argv[1],"r")
genelist=genelist.readlines()

countsubset=0

# while countsubset<500: #500 bootstrap
#     subset = []
#     while len(subset)<len(genelist): #subset is the size of the real list
#         gene=random.choice(genelist)
#         subset.append(gene)
#     sortie = open(sys.argv[2]+"_"+str(countsubset), "a")
#     countsubset += 1
#     for elmt in subset:
#         sortie.write(elmt)


subset = []
sortie = open(sys.argv[2], "w")
while len(subset)<len(genelist): #subset is the size of the real list
    gene=random.choice(genelist)
    subset.append(gene)
for elmt in subset:
    sortie.write(elmt)