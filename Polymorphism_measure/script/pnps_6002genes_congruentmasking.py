from Bio import SeqIO
import os
import sys
import statistics

#path_data="/beegfs/data/mbastian/"
path_data="/home/mbastian/data/"

genes_list_path=path_data+"Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset//145sp/masked/list6002genes"

out_table_path_full=path_data+"Enard_postbusco/pSwithmask/Enard_pS_masked_GQ150QUAL125_fq02to08_6002genes"

genelistfile=open(genes_list_path)
genelist=list(genelistfile)
genelist=[f.replace("\n","") for f in genelist]
genelist=[f.replace("_phylterandbayescode_filtred_145sp.fasta","") for f in genelist]
Efastalist=open(path_data + "Enard_postbusco/1genus_144splist_E","r") #sp
Efastalist=Efastalist.readlines()

sortie_full=open(out_table_path_full,"w")
sortie_full.write("specie\tpS\tpN\tpN/pS\tNSobs\tSobs\tlencallablepos_codon\tnbgenewithSNP\tnbgenetot\n")

lencoding_file=open(path_data+"/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/table_sp2lencallcoding_6002genes","r")
lencoding=lencoding_file.readlines()
sp2lencallable=dict()
for a in lencoding:
    info=a.split()
    call=info[1]
    sp=info[0]
    sp2lencallable[sp]=call


for elmt in Efastalist: #for each sp
    geneprstinlist=0
    genefileexist=0
    seqnotfound = 0
    name = elmt[:-3]
    if os.path.exists(path_data + "VCF_Enard/vcf_annoted/callablesnp/"+name+"/Enard_mam_"+name+"_coding_homofiltred_GQ150QUAL125_fq02to08_6002genes.vcf"):
        #print()
        #print(name,"in progress")
        vcffile = open(path_data + "VCF_Enard/vcf_annoted/callablesnp/"+name+"/Enard_mam_"+name+"_coding_homofiltred_GQ150QUAL125_fq02to08_6002genes.vcf","r")
        line = vcffile.readline()
        gene2NS=dict()
        obsNS=0
        obsS=0
        gene2S=dict()
        gene2length=dict() #used only to count non redudant genes
        genewithsnp=list()
        while line !="":
            info=line.split()
            gene=str(info[4])
            mut=info[7]
            if gene in genelist:
                #if gene not in sp2genefiltred[name]: #so the gene pass the phylter filtre
                #path=str(path_data + "Enard_postbusco/Genes/"+gene+"_no_filter/" + gene+".fasta") #didn't use the length of the gene, just want to know is the seq still here after phylter
                #if os.path.exists(path):
                #    ali = SeqIO.to_dict(SeqIO.parse(path, "fasta")) #really usefull to open the gene ?
                #    if name+"_E" in ali :
                geneprstinlist += 1
                gene2length[gene]=1#used only to count non redudant genes
                if gene not in genewithsnp:
                    genewithsnp.append(gene)
                    genefileexist +=1
                if mut=="NS":
                    if gene in gene2NS.keys():
                        gene2NS[gene]+=1
                    else:
                        gene2NS[gene]=1
                if mut =="S":
                    if gene in gene2S.keys():
                        gene2S[gene]+=1
                    else:
                        gene2S[gene]=1
            line = vcffile.readline()
        nbgene=len(gene2length)
        #print("nb snp in vcf file corresponding to busco genes: ", geneprstinlist, "\nnb busco genes in vcf file:  ",  nbgene)
        if geneprstinlist != 0: #if their is at least one gene with at least one snp
            for gene in genelist:
                if gene not in genewithsnp: #genes from busco without snp in the vcf
                    #if gene not in sp2genefiltred[name]:
                    #if os.path.exists(path_data + "Enard_postbusco/Genes/"+gene+"_no_filter/" + gene+".fasta"):
                    #    ali = SeqIO.to_dict(SeqIO.parse(path_data + "Enard_postbusco/Genes/"+gene+"_no_filter/" + gene+".fasta", "fasta"))
                    #    if name+"_E" in ali:
                            #nbrcodon = int(len(ali[name+"_E"].seq) / 3)
                            #gene2length[gene]=nbrcodon
                        nbgene+=1

            totlen=int(sp2lencallable[name]) /3 #have to be in codon
            for gene, value in gene2NS.items():
                obsNS+=value
            for gene, value in gene2S.items():
                obsS+=value

            PS=float(obsS/totlen)
            PN = float(obsNS / (totlen * 2))  # 2 nuc sur 3 sont NS
            pnps=float(PN/PS)

            sortie_full.write(name+"\t"+str(PS)+"\t"+str(PN)+"\t"+str(pnps)+"\t"+str(obsNS)+"\t"+str(obsS)+"\t"+str(totlen)+"\t"+str(genefileexist)+"\t"+str(nbgene)+"\n")
            print(name+"\t"+str(PS)+"\t"+str(PN)+"\t"+str(pnps)+"\t"+str(obsNS)+"\t"+str(obsS)+"\t"+str(totlen)+"\t"+str(genefileexist)+"\t"+str(nbgene))

        else:
            print(name + " no SNP in busco genes... Have to manually look this vcf")
        vcffile.close()
    else:
        print(name+" no vcf file \n")

    #print("nb of seq not found but gene found in vcf: ", seqnotfound)
