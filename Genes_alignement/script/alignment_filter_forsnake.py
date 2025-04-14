from Bio import SeqIO
import sys
import statistics
import os
import ali

noninformative_characters = "?X-*"
gene=sys.argv[2]
pathgene=sys.argv[1]

#write in fasta format
seqphy=ali.SequenceAlignment(file_name=pathgene,format = "phylip" ) #the seq in output of bmge are in phylip format, have to rewrite in fasta
seqphy.write_fasta_to_file(sys.argv[2], force=True) #write the gene in fasta
del seqphy

listtaxtodel=list()
list_nbinfosite=list()
ali= SeqIO.to_dict((SeqIO.parse(sys.argv[2], "fasta")))

for sp in ali.keys():
    count_noninf=0
    seq=str(ali[sp].seq)
    lenseq = len(seq)
    if lenseq < 10: #justtoremove empty files
        listtaxtodel.append(sp)
    else:
        for nuc in seq:
            if nuc in noninformative_characters:
                count_noninf+=1
        info_site=len(seq)-count_noninf
        if count_noninf!=0:
            frac_noninf = (100 * count_noninf) / lenseq
        else:
            frac_noninf=0
        if frac_noninf > 75:
            listtaxtodel.append(sp)
            #print(gene, " ", sp, " seq removed")
        else:
            list_nbinfosite.append(info_site)
#at this point, have a list of sp to remove from the ali

if listtaxtodel !=0 :
    listtaxtodel=set(listtaxtodel)
    for elmt in listtaxtodel:
        del ali[elmt]
    len_aliafter=len(ali)
    if len(ali)>=146: #80% of 183 as previoulsy used
        sortie = open(sys.argv[3], "w")
        for sp in ali.keys():
            sortie.write(">"+sp+"\n"+str(ali[sp].seq)+"\n")
    else:
        sortie = open(sys.argv[3], "w")
        sortie.write("EMPTY")





