import numpy as np
import statistics

sp2pS=dict()
sp2pN=dict()
sp2pNpS=dict()

splist=open("/home/mbastian/data/Enard_postbusco/1genus_144splist_E","r")
splist=splist.readlines()

#initialise the dict, will contain a list of 500 pS or pN or pN/pS
for elmt in splist:
    elmt=elmt[:-3]
    sp2pS[elmt]=[]
    sp2pN[elmt] = []
    sp2pNpS[elmt] = []

#complete the dico, subset per subset
i=1
while i <500 :
    subset_table=open("/home/mbastian/data/Enard_postbusco/bootstrap_ps/6002genes/subset_"+str(i)+"_pS_6002genes_masked","r")
    line=subset_table.readline() #header
    line = subset_table.readline()
    while line !="":
        (sp, pS, pN, pNpS, length)=line.split()
        if pNpS !="NA":
            sp2pS[sp].append(float(pS))
            sp2pN[sp].append(float(pN))
            sp2pNpS[sp].append(float(pNpS))
        line = subset_table.readline()

    subset_table.close()
    i+=1


#compute the full pS and pNpS
fullpS=open("/home/mbastian/data/Enard_postbusco/bootstrap_ps/6002genes/fulltable_pS_6002genes_masked","r")
sp2fullpS=dict()
sp2fullpnps=dict()
sp2fullpN=dict()
line=fullpS.readline()
line=fullpS.readline()
while line !="":
    (sp, pS, pN, pNpS, length)=line.split()
    if pNpS !="NA":
        sp2fullpS[sp]=float(pS)
        sp2fullpN[sp] = float(pN)
        sp2fullpnps[sp]=float(pNpS)
    line=fullpS.readline()


#compute quantiles 0.25 and 0.75
out=open("/home/mbastian/data/Enard_postbusco/bootstrap_ps/6002genes/bootstrap_quantile_psandpnps_6002genes_masked","w")
out.write("sp\tpS\tq1_pS\tq3_pS\tpN\tq1_pN\tq3_pN\tpNpS\tq1_pNpS\tq3_pN\pS\n")
for sp in splist:
    if sp[:-3] in sp2fullpnps.keys():
        pS_1=np.quantile(sp2pS[sp[:-3]], .25)
        pS_3 = np.quantile(sp2pS[sp[:-3]], .75)

        pN_1=np.quantile(sp2pN[sp[:-3]], .25)
        pN_3 = np.quantile(sp2pN[sp[:-3]], .75)

        pNpS_1=np.quantile(sp2pNpS[sp[:-3]], .25)
        pNpS_3 = np.quantile(sp2pNpS[sp[:-3]], .75)

        #print(sp[:-1]+"\t"+str(sp2fullpS[sp[:-1]])+"\t"+str(pS_1)+"\t"+str(pS_3)+"\t"+str(sp2fullpnps[sp[:-1]])+"\t"+str(pNpS_1)+"\t"+str(pNpS_3) )
        out.write(sp[:-3]+"\t"+str(sp2fullpS[sp[:-3]])+"\t"+str(pS_1)+"\t"+str(pS_3)+"\t"+str(sp2fullpN[sp[:-3]])+"\t"+str(pN_1)+"\t"+str(pN_3)+"\t"+str(sp2fullpnps[sp[:-3]])+"\t"+str(pNpS_1)+"\t"+str(pNpS_3)+"\n")
    else:
        print(sp[:-3])

#determine the 5% outlier data
out=open("/home/mbastian/data/Enard_postbusco/bootstrap_ps/6002genes/bootstrap_q5%_psandpnps_6002genes_masked","w")
out.write("sp\tpS\tq025_pS\tq975_pS\tpN\tq025_pN\tq975_pN\tpNpS\tq025_pNpS\tq975_pNpS\n")
for sp in splist:
    if sp[:-3] in sp2fullpnps.keys():
        pS_1=np.quantile(sp2pS[sp[:-3]], .025)
        pS_3 = np.quantile(sp2pS[sp[:-3]], .975)

        pN_1=np.quantile(sp2pN[sp[:-3]], .025)
        pN_3 = np.quantile(sp2pN[sp[:-3]], .975)

        pNpS_1=np.quantile(sp2pNpS[sp[:-3]], .025)
        pNpS_3 = np.quantile(sp2pNpS[sp[:-3]], .975)

        #print(sp[:-1]+"\t"+str(sp2fullpS[sp[:-1]])+"\t"+str(pS_1)+"\t"+str(pS_3)+"\t"+str(sp2fullpnps[sp[:-1]])+"\t"+str(pNpS_1)+"\t"+str(pNpS_3) )
        out.write(sp[:-3]+"\t"+str(sp2fullpS[sp[:-3]])+"\t"+str(pS_1)+"\t"+str(pS_3)+"\t"+str(sp2fullpN[sp[:-3]])+"\t"+str(pN_1)+"\t"+str(pN_3)+"\t"+str(sp2fullpnps[sp[:-3]])+"\t"+str(pNpS_1)+"\t"+str(pNpS_3)+"\n")
    else:
        print(sp[:-3])
