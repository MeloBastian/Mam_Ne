import os
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
import shutil
bayes=open("/home/mbastian/data/Enard_postbusco/phylter_vs_bayescode/fulldataset/bayescode_outlier_synZsup5","r")
phyl=open("/home/mbastian/data/Enard_postbusco/phylter_vs_bayescode/fulldataset/145sp/phylter_out","r")

b=bayes.readlines()
p=phyl.readlines()
#p=p[4:]

gene2sp=dict() #for each gene, contain the sp to remove
for e in p:
    e=e[:-1]
    data=e.split()
    gene=data[0]
    sp=data[1]
    if data[0] not in gene2sp.keys():
        gene2sp[gene]=[sp]
    else:
        if sp not in gene2sp[gene]:
            gene2sp[gene].append(sp)
        else:
            print("already prst")

for e in b:
    e = e[:-1]
    data = e.split()
    gene=data[0]
    sp=data[1]
    if data[0] not in gene2sp.keys():
        gene2sp[gene]=[sp]
    else:
        if sp not in gene2sp[gene]:
            gene2sp[gene].append(sp)
        else:
            print("already prst")

splist_file = open("/home/mbastian/data/Enard_postbusco/1genus_145splist", "r")
splist = splist_file.readlines()
splist=[f.replace("\n","_E") for f in splist]

genefile=open("/home/mbastian/data/Enard_postbusco/Genes_to_analyse/Genes_filtred/145spgenes_masked/list6007_nonempty", "r")
gene=genefile.readline()[:-1]

countchanges=0
countunchanges=0
countremove=0
while gene !="":
    if os.path.exists("/home/mbastian/data/Enard_postbusco/Genes_to_analyse/Genes_filtred/145spgenes_masked/"+gene):
        geneshort=gene.split("_masked_length")[0]
        if geneshort in gene2sp.keys(): #si il y a des seq à enelver, peut inclure les genes entierement suprimmés par phylter (toutes les seq du genes sont à filtrer) --> donc creer genes vides ?
            new_gene=[] #futur gene
            with open("/home/mbastian/data/Enard_postbusco/Genes_to_analyse/Genes_filtred/145spgenes_masked/"+gene) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id not in gene2sp[geneshort] and record.id in splist: #si la sequence n'est pas à enlever & dans la liste d'especes
                        new_gene.append(record)
            new_alignment = MultipleSeqAlignment(new_gene)
            if len(new_alignment)>116: # 80% of 145
                SeqIO.write(new_alignment, "/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/"+geneshort+"_phylterandbayescode_filtred_145sp.fasta", "fasta")
            countchanges+=1
        else: #si rien à changer, doit qd meme enlever sp pas dans liste
            new_gene = []  # futur gene
            with open("/home/mbastian/data/Enard_postbusco/Genes_to_analyse/Genes_filtred/145spgenes_masked/" + gene) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    if record.id in splist:  # si la sequence est dans la liste d'especes
                        new_gene.append(record)
            new_alignment = MultipleSeqAlignment(new_gene)
            if len(new_alignment) > 116:  # 80% of 145
                SeqIO.write(new_alignment,"/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/" + geneshort + "_phylterandbayescode_filtred_145sp.fasta","fasta")

            #chemin_source = "/home/mbastian/data/Enard_postbusco/Genes_to_analyse/Genes_filtred/"+gene
            #chemin_destination = "/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/"+geneshort+"_phylterandbayescode_filtred_145sp.fasta"

            # Copier le fichier source vers le fichier destination en le renommant
            #shutil.copy(chemin_source, chemin_destination)
            countunchanges+=1
    else:
        print("doesnt exist")

    gene = genefile.readline()[:-1]

# count=0
# for v in gene2sp.values():
#     if len(v)>16:
#         count+=1

# tableout=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/sp2gene_filter", "w")
# for key, value in gene2sp.items():
#     for sp in value:
#         tableout.write(str(sp[:-2])+"\t"+str(key)+"\n") #remove "_E"