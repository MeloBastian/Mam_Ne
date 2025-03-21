from Bio import SeqIO

splist_file = open("/home/mbastian/data/Enard_postbusco/1genus_144splist_E", "r")
splist = splist_file.readlines()
sp2genefilter_file = open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/table_idgenebysp", "r") #remaining genes in each sp, from count_seq.py
sp2genefilter = sp2genefilter_file.readlines()
sp2gene=dict()
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

sp2lencall=dict()

genelist_file=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/list6002genes","r")
genelist=genelist_file.readlines()
genelistshort=[]
for elmt in genelist:
    genelistshort.append(elmt.split("_")[0])

for gene in genelistshort:
    with open("/home/mbastian/data/Enard_postbusco/Genes_to_analyse/"+gene+"binmask.fasta") as mask:
        for record in SeqIO.parse(mask, "fasta"):
            sp =record.id+"_E\n"
            if sp in splist and gene in sp2gene[record.id]:
                seq=list(str(record.seq))
                seqint=[int(element) for element in seq]
                countcall=sum(seqint)
                if record.id not in sp2lencall:
                    sp2lencall[record.id]=int(countcall)
                else:
                    sp2lencall[record.id]+=int(countcall)

out=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/table_sp2lencallcoding_6002genes","w")
for sp, count in sp2lencall.items():
    out.write(sp+"\t"+str(count)+"\n")