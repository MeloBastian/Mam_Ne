from Bio import SeqIO
countseq=0
sp2seq=dict()
a=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/list6002genes","r")
a=a.readlines()
sp2len=dict()

for elmt in a:
    gene=elmt.split("_")[0]
    path="/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/"+elmt[:-1]


    record_dict = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
    for sp in record_dict.keys():
        seq=str(record_dict[sp].seq)
        seq2=seq.replace("-", "")

        if sp not in sp2seq:
            sp2seq[sp]=[elmt[:-1]]
            sp2len[sp]=len(seq2)
        else:
            sp2seq[sp].append(elmt[:-1])
            sp2len[sp] += len(seq2)

        countseq+=1

out=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/table_seqlenbysp", "w")
for key in sp2len.keys():
    out.write(key+"\t"+str(sp2len[key])+"\n")
out2=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/table_idgenebysp", "w")
for key in sp2seq.keys():
    out2.write(key+"\t"+str(sp2seq[key])+"\n")


b=open("/home/mbastian/data/Enard_postbusco/phylter_vs_bayescode/fulldataset/bayescode_outlier_synZsup5", "r")
b=b.readlines()
totgenetoremove=len(b)
for elmt in b:
    gene=elmt.split("\t")[0]+"10andgap75filtred_1genus.fasta"
    sp=elmt.split("\t")[1][:-1]
    if sp not in sp2seq :
        print(sp)
    if sp in sp2seq:
        if gene in sp2seq[sp]:
            totgenetoremove-=1

c=open("/home/mbastian/data/Enard_postbusco/phylter_vs_bayescode/fulldataset/145sp/phylter_out","r")
c=c.readlines()
totgenetoremove=len(c)
for elmt in c:
    gene=elmt.split("\t")[0]
    sp=elmt.split("\t")[1][:-1]
    if sp not in sp2seq :
        print(sp)
    if sp in sp2seq:
        if gene in sp2seq[sp]:
            totgenetoremove-=1