splist_file=open("/home/mbastian/data/Enard_postbusco/1genus_144splist_E","r")
splist=splist_file.readlines()
splist.remove("Homo_sapiens_E\n")
sp2rm=dict()
sp2tot=dict()

for elmt in splist:
    countrm=0
    countall=0
    scaf2pos=dict() #snp to remove
    sp = elmt[:-3]
    print(sp)
    table=open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/"+sp+"/table_coverage_snp_V2_homofiltred_qualiinfo_6002genes", "r")
    line=table.readline()
    line = table.readline()
    line = table.readline()
    while line != "":
        info=line.split()
        fi=int(info[5])/int(info[4])
        if fi >0.8 or fi < 0.2:
            scaf=info[1]
            pos=info[2]

            if scaf not in scaf2pos:
                scaf2pos[scaf]=[pos]
            else:
                scaf2pos[scaf].append(pos)

        line = table.readline()
    vcf = open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/" + sp + "/Enard_mam_"+sp+"_coding_homofiltred_GQ150QUAL125_6002genes.vcf","r")
    out = open("/home/mbastian/data/VCF_Enard/vcf_annoted/callablesnp/" + sp + "/Enard_mam_" + sp + "_coding_homofiltred_GQ150QUAL125_fq02to08_6002genes.vcf","w")
    l=vcf.readline()
    while l != "":
        countall+=1
        info=l.split()
        scaf=info[0]
        pos=info[1]
        if scaf not in scaf2pos:
            out.write(l)
        else:
            if pos not in scaf2pos[scaf]:
                out.write(l)
            else:
                countrm+=1
        l=vcf.readline()
    sp2rm[sp]=countrm
    sp2tot[sp]=countall

sp2raterm=dict()
for sp , rm in sp2rm.items():
    sp2raterm[sp]=(rm/sp2tot[sp])*100

for key, value in sp2raterm.items():
    if value > 10 :
        print(key, value)



