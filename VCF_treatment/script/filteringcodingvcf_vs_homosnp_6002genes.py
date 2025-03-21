splist_file = open("/home/mbastian/data/Enard_postbusco/1genus_144splist_E", "r")
splist = splist_file.readlines()
out_table=open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/table_summary_filtering_snpquali_6002genes","w")
out_table.write("sp\tnbsnptot\tnbsnpquali\trateremovedsnp\trateremovedNS\trateremovedS\n")

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

genelist_file=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/list6002genes","r")
genelist=genelist_file.readlines()
genelistshort=[]
for elmt in genelist:
    genelistshort.append(elmt.split("_")[0])

for elmt in splist:
    countfiltred=0
    count_filtredNS=0
    countNS=0
    count_filtredS=0
    countS=0
    counttot=0
    sp = elmt[:-3]
    print(sp)
    try:
        vcf_coding = open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/" + sp + "/Enard_mam_" + sp + "_coding.vcf", "r")
        line = vcf_coding.readline()
        if line != "":
            out = open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/" + sp + "/Enard_mam_" + sp + "_coding_homofiltred_GQ150QUAL125_6002genes.vcf","w")
            while line != "":
                info = line.split()
                gene=info[4]

                if gene in sp2gene[sp]:
                    G = info[-1].split(":")[0]
                    qual=info[-3]
                    PL_info=info[-1].split(":")[1].split(",")
                    del(PL_info[1])
                    GQ=min(PL_info)
                    type=info[7]
                    if type=="S":
                        countS+=1
                    else:
                        countNS+=1

                    if G != "1/1" and G != "0/0":
                        if int(GQ)>=150 and float(qual)>=125:
                            out.write(line)
                        else:
                            countfiltred+=1
                            if type=="S":
                                count_filtredS+=1
                            else:
                                count_filtredNS+=1
                    counttot+=1
                line = vcf_coding.readline()
            rate_removedsnp=(countfiltred*100)/counttot
            rate_removedNS=(count_filtredNS*100)/countNS
            rate_removedS = (count_filtredS * 100) / countS
            out_table.write(sp+"\t"+str(counttot)+"\t"+str(counttot-countfiltred)+"\t"+str(rate_removedsnp)+"\t"+str(rate_removedNS)+"\t"+str(rate_removedS)+"\n")
    except FileNotFoundError:
        print(sp + " not found")
