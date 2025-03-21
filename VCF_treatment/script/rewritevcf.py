#file=open("/home/mbastian/home/Scripts/Enard_VCF/samplevcf", "r")
#file=open("/home/mbastian/data/VCF_Enard/test/Acomys_cahirinus.raw.vcf")

import sys
file=open("/beegfs/project/nega/heterozygosity/xdisk/denard/denard/mammalian_genomes/vcf_files_lbbe/"+sys.argv[1],"r")
snp=file.readline()
name=sys.argv[1][:-8]
out=open("/beegfs/data/mbastian/VCF_Enard/vcf_D.Enard/"+name+"_2vcftool","a")

while snp !="":
    if snp[0]=="#":
        out.write(str(snp))
    else:
        snp=snp.split()
        if len(snp)==10:
            score=snp[7].split(";")
            DP=score[0].split("=")
            if DP[0]=="DP":
                snp[-2]=snp[-2]+":"+DP[0]
                snp[-1]=snp[-1]+":"+DP[1]
                out.write(str(snp[0])+"\t"+str(snp[1])+"\t"+str(snp[2])+"\t"+str(snp[3])+"\t"+str(snp[4])+"\t"+str(snp[5])+"\t"+str(snp[6])+"\t"+str(snp[7])+"\t"+str(snp[8])+"\t"+str(snp[9])+"\n")

        else:
            print(snp)
    snp = file.readline()
print("done" + name)