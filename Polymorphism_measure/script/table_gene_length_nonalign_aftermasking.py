import os
import pandas as pd
from Bio import SeqIO

#path_data="/beegfs/data/mbastian/"
path_data="/home/mbastian/data/"

#gene list, not to search gene in data but juste to have a list
genes_list_path=path_data+"Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/list6002genes"
genelistfile=open(genes_list_path)
genelist=genelistfile.readlines()
genelist=[f.replace("\n","") for f in genelist]
genelist=[f.replace("_phylterandbayescode_filtred_145sp.fasta","") for f in genelist] #only to have a genelist, indep of aligned or not

#sp list
Efastalist=open(path_data + "Enard_postbusco/1genus_144splist_E","r")
Efastalist=Efastalist.readlines()
splist=[]
for e in Efastalist:
    splist.append(e[:-3])
splist.remove("Homo_sapiens")

#removedseqbyfiltering
#removedseq_file=open(path_data+"Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/sp2gene_filter", "r") #sp\tgeneid
#removedseq=removedseq_file.readlines()
#sp2removedseq=dict()
#for elmt in removedseq:
#    info=elmt.split()
#    if info[0] not in sp2removedseq.keys():
#        sp2removedseq[info[0]]=[info[1]]
#    else:
#        sp2removedseq[info[0]].append((info[1]))


genelist_file = open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/list6002genes","r")
genelist = genelist_file.readlines()
genelistshort = []
for elmt in genelist:
    genelistshort.append(elmt.split("_")[0])

sp2genefilter_file = open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/table_idgenebysp", "r") #remaining genes in each sp, from count_seq.py
sp2genefilter = sp2genefilter_file.readlines()
sp2gene= {}
for line in sp2genefilter:
    info=line.split()
    if info[0] not in sp2gene.keys():
        sp2gene[info[0][:-2]]=[]
    for elmt in info[1:]:
        elmt=elmt.replace("[", "")
        elmt=elmt.replace("]", "")
        elmt = elmt.replace(",", "")
        elmt=elmt.replace("_phylterandbayescode_filtred_145sp.fasta'", "")
        gene=elmt.replace("'","")
        sp2gene[info[0][:-2]].append(gene)

#have to compute len callable position of genes by sp --> to count nb of 1 in binmaskedgenes
gene2sp2length= {}
#sp=splist[0]
for gene in genelistshort:
    sp2length={}
    with open("/home/mbastian/data/Enard_postbusco/Genes_to_analyse/"+gene+"binmask.fasta") as mask:
        for record in SeqIO.parse(mask, "fasta"):
            sp =record.id+"_E"
            if record.id in splist and gene in sp2gene[record.id]: #the sq for this espece is selected
                seq=list(str(record.seq))
                seqint=[int(element) for element in seq]
                countcall=sum(seqint)
                sp2length[record.id]=int(countcall)

    gene2sp2length[gene]=sp2length

sp2gene2length = {}
for gene, especes in gene2sp2length.items():
    # Parcourir les sous-dictionnaires espece:valeur
    for espece, valeur in especes.items():
        # Vérifier si l'espece existe déjà dans le dictionnaire reformatté
        if espece not in sp2gene2length:
            # Si non, ajouter une nouvelle entrée avec un nouveau dictionnaire
            sp2gene2length[espece] = {}

        # Ajouter une entrée (gene, valeur) au dictionnaire correspondant à l'espece
        sp2gene2length[espece][gene] = valeur


#write spxgene table for gene length
data=[]
for sp in splist:
    line=[]
    # if sp in sp2removedseq.keys():
    #     for gene in genelist:
    #         if gene in sp2removedseq[sp]: #si le gene a été filtré pour cette esp
    #             line.append("NA")
    #         else:
    #             try:
    #                 line.append(sp2gene2length[sp][gene])
    #             except:
    #                 line.append(0) #si l'espece ne possede pas ce gene
    # else:
    for gene in genelistshort:
        try:
            line.append(sp2gene2length[sp][gene])
        except:
            line.append(0) #si l'espece ne possede pas ce gene

    data.append(line)
df=pd.DataFrame(data, columns=genelistshort)
df.insert(loc=0, column="species", value=splist)

df.to_csv(path_data+"Enard_postbusco/coding_callable/withgeneinfo/table_spxgenecallablelength_6002genes_fq02to08_masked")


#compute S and NS by gene
sp2gene2S={}
sp2gene2NS={}
for sp in splist:
    vcffile = open(path_data + "VCF_Enard/vcf_annoted/callablesnp/" + sp + "/Enard_mam_" + sp + "_coding_homofiltred_GQ150QUAL125_fq02to08_6002genes.vcf", "r")
    print(sp, "in progress")
    line = vcffile.readline()
    gene2NS = {}
    obsNS = 0
    obsS = 0
    gene2S = {}
    genewithsnp = []
    while line != "":
        info = line.split()
        gene = str(info[4])
        mut = info[7]
        if gene in genelistshort:
            if gene not in genewithsnp:
                genewithsnp.append(gene)
            if mut == "NS":
                if gene in gene2NS.keys():
                    gene2NS[gene] += 1
                else:
                    gene2NS[gene] = 1
            if mut == "S":
                if gene in gene2S.keys():
                    gene2S[gene] += 1
                else:
                    gene2S[gene] = 1
        line = vcffile.readline()
    sp2gene2NS[sp]=gene2NS
    sp2gene2S[sp]=gene2S

#write S and NS table
dataS=[]
dataNS=[]
for sp in splist:
    lineS=[]
    lineNS=[]
    # if sp in sp2removedseq.keys():
    #     for gene in genelist:
    #         if gene not in sp2removedseq[sp]:
    #             try:
    #                 lineS.append(sp2gene2S[sp][gene])
    #             except:
    #                 lineS.append(0) #0 bc possible to have no snp in a gene, not a NA
    #             try:
    #                 lineNS.append(sp2gene2NS[sp][gene])
    #             except:
    #                 lineNS.append(0)
    #         else:
    #             lineS.append("NA") #bc gene filtred
    #             lineNS.append("NA")
    #else:
    for gene in genelistshort:
        try:
            lineS.append(sp2gene2S[sp][gene])
        except:
            lineS.append(0)
        try:
            lineNS.append(sp2gene2NS[sp][gene])
        except:
            lineNS.append(0)
    dataS.append(lineS)
    dataNS.append(lineNS)
dfS=pd.DataFrame(dataS, columns=genelistshort)
dfS.insert(loc=0, column="species", value=splist)
dfNS=pd.DataFrame(dataNS, columns=genelistshort)
dfNS.insert(loc=0, column="species", value=splist)

dfS.to_csv(path_data+"Enard_postbusco/bootstrap_ps/6002genes/table_genexsp_Scount_6002genes_fq02to08_masked")
dfNS.to_csv(path_data+"Enard_postbusco/bootstrap_ps/6002genes/table_genexsp_NScount_6002genes_fq02to08_masked")



#to do : remove filtred seq by adding NA : create a dico sp2removedseq and verify to not be in this dico

