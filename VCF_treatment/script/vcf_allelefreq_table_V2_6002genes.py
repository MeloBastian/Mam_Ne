import os
import random
import Bio

genelist_file=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/list6002genes","r")
genelist=genelist_file.readlines()
genelist = [f.replace("_phylterandbayescode_filtred_145sp.fasta\n", "") for f in genelist]
sp2genefilter_file = open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/table_idgenebysp", "r") #remaining genes
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



splist_file=open("/home/mbastian/data/Enard_postbusco/1genus_144splist_E","r")
splist=splist_file.readlines()
splist.remove("Homo_sapiens_E\n")

listrawvcf = os.listdir("/home/mbastian/data/VCF_Enard/vcf_D.Enard/")
listrawvcf_2vcftool=[]
for elmt in listrawvcf:
    if "2vcftool"  in elmt:
        listrawvcf_2vcftool.append(elmt)

for elmt in splist:
    sp=elmt[:-3]
    print(sp)

    try:
        coding_vcf_file=open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/"+sp+"/Enard_mam_"+sp+"_coding_homofiltred_GQ150QUAL125_6002genes.vcf","r")
    except FileNotFoundError:
        print(sp+ " not found")

    for elmt in listrawvcf_2vcftool:
        if sp in elmt:
            raw_vcf_file = open("/home/mbastian/data/VCF_Enard/vcf_D.Enard/"+elmt, "r")

    #recup info dp4 dans un dict
    scafpos2dp4=dict()
    raw_line=raw_vcf_file.readline()
    while raw_line != "":
        if raw_line[0]!="#": #NW_020835294.1	1979	.	C	T	222.239	.	DP=87;VDB=0.506048;SGB=-0.693136;RPBZ=0.811793;BQBZ=-1.38749;SCBZ=-0.136852;FS=0;MQ0F=0;AC=1;AN=2;DP4=23,25,17,18;MQ=60	GT:PL:DP	0/1:255,0,255:87
            info=raw_line.split()
            scafpos=info[0]+"\t"+info[1]
            dp4=info[7].split(";")[-2].split("=")[1]
            scafpos2dp4[scafpos]=dp4
        raw_line = raw_vcf_file.readline()
    raw_vcf_file.close()


    out = open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/"+sp+"/table_coverage_snp_V2_homofiltred_qualiinfo_6002genes", "w")
    #out = open("/home/mbastian/Sauv-sept2023/" + sp + "/table_coverage_snp_homofiltred", "w")
    out.write("#use DP4 : nb of high quality ref-forward, ref-reverse, alt-forward and alt-reverse. Sum DP4=mi<DP because use only 'high quality'. ki=sum ref \n")
    out.write("gene\tscaffold\tposition\ttype\tmi\tki_perm\tqt_ref\tqt_ref_forward\tqt_ref_reverse\tqt_alt\tqt_alt_forward\tqt_alt_reverse\tDP\tGQ\tqual\n")
    coding_line = coding_vcf_file.readline()
    while coding_line != "":
        info=coding_line.split()
        gene = info[4]
        if gene in sp2gene[sp]: #some sp are removed from some genes by phylter analyse
            scaf=info[0]
            pos=info[1]
            type=info[7]
            qual=info[-3]
            PL_info=info[-1].split(":")[1].split(",")
            del(PL_info[1])
            GQ=min(PL_info)

            dp4=scafpos2dp4[scaf+"\t"+pos]
            cov=dp4.split(",")
            #don't know how ref is choosed so random permutation for ki or mi-ki=li
            permutation = bool(random.choice([True, False]))

            ki=int(cov[0])+int(cov[1]) #choose if true
            mi=int(cov[0])+int(cov[1])+int(cov[2])+int(cov[3])
            li=mi-ki #choose if false

            DP=info[-1].split(":")[-1]
            if permutation == True:
                out.write(gene+"\t"+scaf+"\t"+pos+"\t"+type+"\t"+str(mi)+"\t"+str(ki)+"\t"+str(ki)+"\t"+cov[0]+"\t"+cov[1]+"\t"+str(li)+"\t"+cov[2]+"\t"+cov[3]+"\t"+str(DP)+"\t"+str(GQ)+"\t"+str(qual)+"\n")
            else:
                out.write(gene+"\t"+scaf+"\t"+pos+"\t"+type+"\t"+str(mi)+"\t"+str(li)+"\t"+str(ki)+"\t"+cov[0]+"\t"+cov[1]+"\t"+str(li)+"\t"+cov[2]+"\t"+cov[3]+"\t"+str(DP)+"\t"+str(GQ)+"\t"+str(qual)+"\n")
            coding_line = coding_vcf_file.readline()
        else:
            coding_line = coding_vcf_file.readline()
    coding_vcf_file.close()


