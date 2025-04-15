splist_file = open("/home/mbastian/data/Enard_postbusco/1genus_145splist", "r")
splist = splist_file.readlines()
#splist=splist[0:3]
genelist=[]

sp2gene2seq=dict()


print("write masked gene position")
for s in splist:
    sp=s[:-1]
    print(sp)
    #all exons positions
    allexons=open("/home/mbastian/data/Enard_postbusco/Refexons/"+sp+"_E_phase1.gff_sorted","r").readlines()
    gene2exons_scafcoord=dict() #list of tuple containing start and end position of each exon for a given gene, in phase 1
    gene2brin = dict()
    for elmt in allexons: # an elmt is an exon
        info=elmt.split()
        start=int(info[2])
        end=int(info[3])
        gene=info[-1]
        if gene not in genelist: #for the gene writing part
            genelist.append(gene)
        pos=(int(start), int(end))
        if gene not in gene2exons_scafcoord.keys():
            brin = info[4]
            gene2brin[gene]=brin
            gene2exons_scafcoord[gene]=[]
        gene2exons_scafcoord[gene].append(pos)

        if end < start : # always false
            print(info, start, end)


    gene2exon_genecoord=dict() #in phase 1
    for gene, scafcoords in gene2exons_scafcoord.items():
        nbexons=len(scafcoords)

        #exon1, gene start , pos start= 1
        startCDS=scafcoords[0][0]
        endCDS = scafcoords[0][1]
        startgene=1
        endgene=endCDS-startCDS +1
        pos=(startgene, endgene)
        gene2exon_genecoord[gene]=[pos]

        if nbexons > 1:  #si multi-exon
            i=2
            while i < nbexons +1:
                startCDS = scafcoords[i-1][0]
                endCDS = scafcoords[i-1][1]
                startgene=gene2exon_genecoord[gene][-1][-1] +1#end position of the last record exon
                endgene= gene2exon_genecoord[gene][-1][-1] +1 + endCDS - startCDS
                pos=(startgene, endgene)
                gene2exon_genecoord[gene].append(pos)

                if endgene<startgene: #alwawys false
                    print(scafcoords[i-1], pos)
                i+=1


    maskedexon=open("/home/mbastian/data/Enard_postbusco/coding_callable/withgeneinfo/"+sp+"_masked_exon_geneinfo.bed").readlines() #in phase 0
    gene2mask_scafcoord=dict() #phase 0 (bedtools)
    for elmt in maskedexon: # an elmt is an exon
        info=elmt.split()
        start=int(info[1])#+1 #phase0-->phase1
        end=int(info[2])#+1#phase0-->phase1
        gene=info[-1]
        pos=(int(start), int(end))
        if gene not in gene2mask_scafcoord.keys():
            gene2mask_scafcoord[gene]=[]
        gene2mask_scafcoord[gene].append(pos)


    gene2mask_genecoord=dict() #phase 1
    for gene, scafcoords in gene2mask_scafcoord.items():
        gene2mask_genecoord[gene]=[]
        nbexons=len(scafcoords)
        i=1
        while i<nbexons +1:
            SMS=scafcoords[i-1][0] #start mask scaf
            EMS = scafcoords[i-1][1]
            exonnb=0 #codingexon
            while exonnb < len(gene2exons_scafcoord[gene]):
                exon= gene2exons_scafcoord[gene][exonnb]
                SES=exon[0] #start exon scaf
                EES=exon[1]
                if SMS >= SES and EMS>= SES and SMS <= EES and EMS <= EES:# if start and end of masked part is in this exon
                    SEG = gene2exon_genecoord[gene][exonnb][0]
                    EEG = gene2exon_genecoord[gene][exonnb][1]
                    # if gene2brin[gene]=="+":
                    SMG=SEG+(SMS - SES)  # start position of masked exon = start of mask in CDS - start of exon in CDS
                    EMG=SMG+(EMS-SMS)#+1 #end position of masked exon = end of exon in CDS - end of mask in CDS --> hypothesis masked position are nested in the same exon. (true normaly...)
                    pos=(SMG, EMG)
                    # else:
                    #     SMG=SEG+abs(EMS-EES)
                    #     EMG=SMG+abs(SMS-EMS)#+1
                    #     pos=(SMG, EMG)

                    gene2mask_genecoord[gene].append(pos)
                exonnb+=1
            i+=1

    #ecrire une sequence de la taille du CDS, si position pas dans masque, ecrire 1, sinon ecrire 0
    #write a sequence with len equal to the CDS len, if the position is not in the mask, write a 1, else write a 0 --> allow later to sum the "1" and obtain the callable fraction of each CDS
    gene2len=dict()
    for gene in gene2exon_genecoord:
        totlen=0
        for elmt in gene2exon_genecoord[gene]:
            start=elmt[0]
            end=elmt[1]
            length=end-start+1
            totlen+=length
        gene2len[gene]=totlen

    gene2seqbin=dict()
    for gene in gene2len:
        totlen=gene2len[gene]
        a=0
        seq=[]
        while a <=totlen-1:#because gene pos in base 1 and python in base 0
            done = False
            if gene in gene2mask_genecoord:
                for elmt in gene2mask_genecoord[gene]:
                    if a >= elmt[0]-1 and a <=elmt[1]-1 and done == False: 
                        seq.append("0")
                        done=True
            if done ==False:
                seq.append("1")
            a+=1
         if gene2brin[gene]=="-":
             seq="".join(reversed(seq))
         else:
             seq="".join(seq)


        gene2seqbin[gene]=seq
    sp2gene2seq[sp]=gene2seqbin
    #print(seq)

print("write gene binmask of "+str(len(genelist))+" genes")
for gene in genelist:
    out=open("/home/mbastian/data/Enard_postbusco/Genes_to_analyse/"+gene+"binmask.fasta","w")
    sp2seq=dict()
    for sp in sp2gene2seq:
        if gene in sp2gene2seq[sp]:
            out.write(">"+sp+"\n"+sp2gene2seq[sp][gene]+"\n")
            #sp2seq[sp]=sp2gene2seq[sp][gene]
