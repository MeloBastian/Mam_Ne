import sys
import pandas as pd
import numpy as np

#import genexsp table
dflength= pd.read_csv(sys.argv[3])
#dflength=pd.read_csv("/home/mbastian/data/Enard_postbusco/coding_callable/withgeneinfo/table_spxgenecallablelength_6002genes_fq02to08_masked")
dflength.set_index('species', inplace = True)
dfS= pd.read_csv(sys.argv[4])
#dfS=pd.read_csv("/home/mbastian/data/Enard_postbusco/bootstrap_ps/6002genes/table_genexsp_Scount_6002genes_fq02to08_masked")
dfS.set_index('species', inplace = True)
dfNS= pd.read_csv(sys.argv[5])
#dfNS=pd.read_csv("/home/mbastian/data/Enard_postbusco/bootstrap_ps/6002genes/table_genexsp_NScount_6002genes_fq02to08_masked")
dfNS.set_index('species', inplace = True)

#import the subset
subset_file=open(sys.argv[1],"r")
#subset_file=open("/home/mbastian/data/Enard_postbusco/Genes_after_phylterandbayescode_filtering/fulldataset/145sp/masked/list6002genes", "r")
subset=subset_file.readlines()
subset=[f.replace("_phylterandbayescode_filtred_145sp.fasta\n","") for f in subset]


#some genes in subset ar not in lenght, pS and pN table, dont know why
for gene in subset:
    if gene not in dflength or gene not in dfNS or gene not in dfS:
        subset.remove(gene)

#dont know why but some genes didnt remove at first occurence
for gene in subset:
    if gene not in dflength or gene not in dfNS or gene not in dfS:
        subset.remove(gene)

#reduce tables to subsetlist
subsetuniq=list(set(subset))#if use bootstrap subset, redudunt columns difficults to deal with later so we reduce the table
subdflength=dflength[subsetuniq]
subdfS=dfS[subsetuniq]
subdfNS=dfNS[subsetuniq]


out=open(sys.argv[2],"w")
#out=open("/home/mbastian/data/Enard_postbusco/bootstrap_ps/6002genes/fulltable_pS_6002genes_masked", "w")
out.write("specie\tpS\tpN\tpN/pS\tlength_aa\n")



for sp in list(subdflength.index):
    sumNS = 0
    sumS = 0
    sumlengthaa = 0
    for gene in subset: #take count of all genes and their redundance
        #print(gene)
        length=str(subdflength[gene][sp])
        S=str(subdfS[gene][sp])
        NS = str(subdfNS[gene][sp])
        if length!="nan" and NS!="nan" and S!="nan": #remove nan (= filtred seq)
            if length != "0" :#gene exist
                #print(length, S, NS)
                sumNS+=float(NS)
                sumS+=float(S)
                sumlengthaa+=(float(length)/3)
    pS=sumS/sumlengthaa
    pN=sumNS/(sumlengthaa*2)
    try:
        pnps=pN/pS #homosapiens didnt work
    except:
        pnps="NA"
    out.write(sp+"\t"+str(pS)+"\t"+str(pN)+"\t"+str(pnps)+"\t"+str(sumlengthaa)+"\n")
