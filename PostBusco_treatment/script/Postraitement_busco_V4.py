
import os
from Bio import SeqIO
from Bio.Seq import Seq
from ali import translate_nucseq
import sys

split = sys.argv[1]
print("split: " + split)
PathBanque = sys.argv[2]
PathData = sys.argv[3]

especes_list_path = sys.argv[4]
text_file = open(especes_list_path, "r")
lines = text_file.readlines()
lines = [f.replace(".fna\n", "") for f in lines]
especes_list = [f.replace("mam_", "") for f in lines]
gene_list = list()

print(especes_list)

sp2gene2seq = dict()  # each species key contain a gene2sec dict 
sp2gene2ref = dict()  # each species key contain a gene2ref dict
sp2gene2scafold = dict() # each species key contain a gene2scaffold dict
sp2gene2scafoldpos = dict() # each species key contain a gene2scafoldpos dict

##########################################################################################################################
for sp in especes_list:  # For each sp, construct a gene2seq and gene2ref dict wich respectively contain every gene nuc sequences and exons positions, after filtering steps.
    print()
    print(sp)
    genome_dict = SeqIO.to_dict(SeqIO.parse("/beegfs/banque/mbastian/Enard/genomes/"+"mam_" + sp + ".fna", "fasta")) #genome
    path_sp = PathBanque + "/Busco_splited_species/mam_" + sp + ".fna/" 
    gene2seq = dict()  # contain every nuc sequences after filtering
    gene2ref = dict()  # contain exons informations for each genes after filtering
    gene2count = dict()  # count the occurence of each gene. They are all suposed to be seen only once in the total genome.
    gene2aaseq = dict()  # countain the aa sequences 
    listegenesinglesp = list()  # a list of the real single copy genes in the genome
    gene2scafold = dict() #relation between the gene and the scaffold id
    gene2scafoldpos = dict() # relation between the gene and is position in the scaffold

    #some monitoring: 
    compte_seqdiffrerun = 0
    compte_seqdiffinitial = 0 #count the number of genes presenting a nuc sequence with traduction in aa is not comaptible with the aa seq provided
    compte_taille = 0 #count the number of CDS whith a size not multiple of 3
    compte_falsesingle = 0 #count the number of genes pretended single by BUSCO which are not
    compte_exon=0
    compte_exons_out_false=0
    compte_exons_in_false=0
    nocorresp=0
    nbrseq=0
    nbrwrongseq_in=0
    nbrwrongseq_out=0

#a first reading of the differents splits to fill the listegenesinglesp and gene2aaseq
    #singleseqpath="/home/mbastian/Bureau/busco_mam_Hystrix_cristata.fna_0/run_mammalia0_odb10/busco_sequences/single_copy_busco_sequences/"
    singleseqpath = path_sp + "busco_mam_" + sp + ".fna_" + str(split) + "/run_mammalia" + str(split) + "_odb10/busco_sequences/single_copy_busco_sequences/"   #contain a file per gene, flaged as single copy by BUSCO
    for gene in os.listdir(singleseqpath):
        id = gene[:-4]  # rm the ".faa"
        if id not in gene2count:  # sometimes their is multiple occurence of the aa sequences in the file but they are identical. Here its a shortcut to using only the first occurence of the aa seq
            gene2count[id] = 1
            aapath = singleseqpath + id + ".faa"  
            if os.path.isfile(aapath): 
                file = open(aapath, "r")
                record = SeqIO.read(file, "fasta")
                idseq =str(record.id)
                idseq =idseq.split(":")
                idscafold=idseq[0] #scaffold name
                idscafold_pos=idseq[1] #position in the scaffold
                gene2scafold[id]=idscafold
                gene2scafoldpos[id]=idscafold_pos
                gene2aaseq[id] = str(record.seq)
        else:
            gene2count[id] += 1
            if gene2count[id] >1 :
                print("gene" +id+" not unique")
    for keys in gene2count.keys():  # list genes that have only been observed once: count value =1
        if gene2count[keys] == 1:
            listegenesinglesp.append(keys)

# a second reading of the species splits to extract the codon sequences from the genes in the listgenesinglesp and performing different control (the genes id and the sequences are not in the same files)
# Busco seems to run twice as their is an "initial" and a "rerun" dir. Read the two one after the other.
    singleseqpath = path_sp + str(split) + "/run_mammalia_odb10/busco_sequences/single_copy_busco_sequences/"
    fichier = "mam_" + sp + ".fna" + ".codon.fas"  
    #pathinitial="/home/mbastian/Bureau/busco_mam_Hystrix_cristata.fna_0/run_mammalia0_odb10/metaeuk_output/initial_results/mam_Hystrix_cristata.fna.codon.fas"
    pathinitial = path_sp + "busco_mam_" + sp + ".fna_" + str(split) + "/run_mammalia" + str(split) + "_odb10/metaeuk_output/initial_results/" + fichier
    pathrerun = path_sp + "busco_mam_" + sp + ".fna_" + str(split) + "/run_mammalia" + str(split) + "_odb10/metaeuk_output/rerun_results/" + fichier

    #  "initialresults" files
    alreadytreated=list()
    if os.path.exists(pathinitial):  
        initialcodon = open(pathinitial, "r")
        while 1:
            idfull = initialcodon.readline()  # first line = exon info
            if idfull == "":  
                break
            seq = initialcodon.readline().strip()  # second line = sequence
            idlist = idfull.split("_")
            id = idlist[0][1:] 
            id2 = id.split("|") 
            id3 = id2[0]  
            reffull = idfull.split("|", maxsplit=5)
            scafinitial=reffull[1] #scaf ref from initial files
            ref = reffull[-1]  # exon info
            #print(ref)


            if id3 in listegenesinglesp and id3 not in alreadytreated :  # if the gene is in the single copy gene list and the right seq hasn't already been found
                nbrseq += 1
                # the CDS have to be multiple of 3
                if len(seq) % 3 != 0:
                    print("the sequence is not a multiple of 3", sp, id3)
                    compte_taille += 1
                else:
                    seq_up = seq.upper()
                    # the traduction of the codon sequence  in "initial results" have to be the same than the aa seq 
                    translate_nuc = translate_nucseq(seq_up)
                    translate_nuc_up = translate_nuc.upper()
                    aaseq = gene2aaseq[id3]
                    aaseq_up = aaseq.upper()
                    if translate_nuc_up != aaseq_up:
                        compte_seqdiffinitial += 1

                    # compare the exons seq obtained after BUSCO analyse and the genome sequence
                    scafold = gene2scafold[id3]
                    if scafold == scafinitial : #test if the scafold id from busco seq are the same than in the initial file
                        refsplit = ref.split("|")
                        exons = refsplit[3:]  # info position exons
                        seqgenome = genome_dict[scafold].upper()  # the scafold seq
                        errorinfullseq_in = 0
                        # errorinfullseq_out=0
                        CDSfromgenome = ""  # to concatenate the exons froms genome and recreate the CDS in order to compare it with the CDS from Busco
                        for elmt in exons:
                            compte_exon += 1
                            exon = elmt.split(":")
                            start_in = exon[0].split("[")[1][:-1]
                            start_in = int(start_in)
                            end_in = exon[1].split("[")[1][:-1]
                            end_in = int(end_in)   # +1 to have the last position include in the interval in python

                            # keep the exon seq from the scafold with positions
                            if int(start_in) < int(end_in):
                                start = int(start_in)
                                end = int(end_in) + 1  # +1 to have the last position include in the interval in python
                                seq_exon_in = seqgenome[start:end]
                                seq_exon = str(seq_exon_in.seq)

                            else:
                                if int(start_in)> int(end_in):
                                    start = int(end_in)
                                    end = int(start_in) + 1
                                    seq_exon_in = seqgenome[start:end]
                                    seq_exon_in=str(seq_exon_in.seq)
                                    seq_exon_in=Seq(seq_exon_in)
                                    seq_exon_in_reversecomp=seq_exon_in.reverse_complement()
                                    seq_exon = str(seq_exon_in_reversecomp)

                                else:
                                    print("error in position")

                            if len(seq_exon) != 0:  # some exons look like having a size 0
                                if str(seq_exon).upper() not in str(seq_up):
                                    compte_exons_in_false += 1
                                    errorinfullseq_in += 1
                                    nocorresp += 1
                                    print("genome seq_exon not in busco CDS ")
                                else:
                                    CDSfromgenome += str(seq_exon)  # concatenate the exons
                            else:
                                print("len exon = 0") #be careful to the previous +1
                                compte_exons_in_false += 1
                                errorinfullseq_in += 1

                        if errorinfullseq_in != 0:  # after the writing of the CDS, if one exon is not good
                            nbrwrongseq_in += 1
                        else:
                            if CDSfromgenome != seq_up:  # try if the CDS reconstructed from the genome is the same than the one from busco 
                                print("CDS from genome different from Busco")
                                # nbrwrongCDS+=1
                            else:
                                gene2seq[id3] = seq  # If every filters are passed, write the sequence in gene2seq dict
                                gene2ref[id3] = ref
                                alreadytreated.append(id3)  # if a good CDS is found with all the filter corect, do not try to find an other CDS for this gene (often usefull to optimise the "rerun files" reading)
                    # else:
                    #     print("id scafold from initial files is different than busco sequence file ")

    else:
        print(pathinitial+ "doesn't exist")  
#nbrexonfalseout=(100*compte_exons_out_false)/compte_exon
# nbrexonfalsein=(100*compte_exons_in_false)/compte_exon
# tauxwronseq_in=(100*nbrwrongseq_in)/nbrseq
# tauxwronseq_out=(100*nbrwrongseq_out)/nbrseq
# print("initial result : nbr exon false out "+str(nbrexonfalseout))
# print("initial result : nbr exon false in "+str(nbrexonfalsein))
# print("initialresult : nbr sequence false in "+str(tauxwronseq_in))
# print("initialresult : nbr sequence false out "+str(tauxwronseq_out))
#
# compte_falsesingle = 0
# compte_exon = 0
# compte_exons_out_false = 0
# compte_exons_in_false = 0
# nbrseq = 0
# nbrwrongseq_in = 0
# nbrwrongseq_out = 0

    # "rerunresults" files reading, same process as for the "initial results"
    if os.path.exists(pathrerun):  
        initialcodon = open(pathrerun, "r")
        while 1:
            idfull = initialcodon.readline()  
            if idfull == "":  
                break
            seq = initialcodon.readline().strip()  
            idlist = idfull.split("_")
            id = idlist[0][1:]  
            id2 = id.split("|")  
            id3 = id2[0]  
            reffull = idfull.split("|", maxsplit=5)
            ref = reffull[-1]
            #print(ref)


            if id3 in listegenesinglesp and id3 not in alreadytreated : 
                nbrseq += 1

                # the CDS have to be multiple of 3
                if len(seq) % 3 != 0:
                    print("the sequence is not a multiple of 3", sp, id3)
                    compte_taille += 1
                else:
                    seq_up = seq.upper()
                    translate_nuc = translate_nucseq(seq_up)
                    translate_nuc_up = translate_nuc.upper()
                    aaseq = gene2aaseq[id3]
                    aaseq_up = aaseq.upper()
                    if translate_nuc_up != aaseq_up:
                        # print("sequence nuc", sp, id3, "differentes de la sequence aa")
                        compte_seqdiffinitial += 1

                    scafold = gene2scafold[id3]
                    if scafold == scafinitial : 
                        refsplit = ref.split("|")
                        exons = refsplit[3:]  
                        seqgenome = genome_dict[scafold].upper()  
                        errorinfullseq_in = 0
                        # errorinfullseq_out=0
                        CDSfromgenome = ""  
                        for elmt in exons:
                            compte_exon += 1
                            exon = elmt.split(":")
                            start_in = exon[0].split("[")[1][:-1]
                            start_in = int(start_in)
                            end_in = exon[1].split("[")[1][:-1]
                            end_in = int(end_in)

                            if int(start_in) < int(end_in):
                                start = int(start_in)
                                end = int(end_in) + 1  # +1 to have the last position include in the interval in python
                                seq_exon_in = seqgenome[start:end]
                                seq_exon = str(seq_exon_in.seq)

                            else:
                                if int(start_in)> int(end_in):
                                    start = int(end_in)
                                    end = int(start_in) + 1
                                    seq_exon_in = seqgenome[start:end]
                                    seq_exon_in=str(seq_exon_in.seq)
                                    seq_exon_in=Seq(seq_exon_in)
                                    seq_exon_in_reversecomp=seq_exon_in.reverse_complement()
                                    seq_exon = str(seq_exon_in_reversecomp)
                                else:
                                    print("error in position")

                            if len(seq_exon)!=0: #some exons look like having a size 0
                                if str(seq_exon).upper() not in str(seq_up):
                                    compte_exons_in_false += 1
                                    errorinfullseq_in += 1
                                    nocorresp+=1
                                    print("genome seq_exon not in busco CDS ")

                                else:
                                    CDSfromgenome += str(seq_exon)  
                            else:
                                print("len exon = 0") #be careful to the previous +1
                                compte_exons_in_false += 1
                                errorinfullseq_in += 1

                        if errorinfullseq_in != 0:  # after the rewrite of the CDS, if one exon is not well
                            nbrwrongseq_in += 1
                        else:
                            if CDSfromgenome != seq_up:  
                                print("CDS from genome different from Busco")
                                # nbrwrongCDS+=1
                            else:
                                gene2seq[id3] = seq  
                                gene2ref[id3] = ref
                                alreadytreated.append(id3) 
                    # else:
                    #     print("id scafold from rerun files is different than busco sequence file ")

    else:
        print(pathinitial+ "doesn't exist") 
    nbrexonfalsein=(100*compte_exons_in_false)/compte_exon
    tauxwronseq_in=(100*nbrwrongseq_in)/nbrseq
    print("initial+rerun result : % exon false "+str(nbrexonfalsein))
    print("initial+rerunn result : % CDS false  "+str(tauxwronseq_in))
    print("no coresp exon from genome and CDS from busco", str(nocorresp))



    print(sp, "nb of single genes before filtering =", len(listegenesinglesp),", nb of seq kept =", len(gene2seq), " --> rate of kept genes =",(len(gene2seq) / len(listegenesinglesp)) * 100, "%")

    sp2gene2seq[sp] = gene2seq
    sp2gene2ref[sp] = gene2ref
    sp2gene2scafold[sp]=gene2scafold
    sp2gene2scafoldpos[sp]=gene2scafoldpos

    gene_list += listegenesinglesp

gene_list = list(set(gene_list))
sortie = open(PathData + "Genes/gene_list{}".format(str(split)), "w")
sortie.write(str(gene_list))

##########################################################################################################################
compte_smaller = 0

#write fasta files and summary of the list of genes per species and list of specie per gene.
for gene in gene_list:  # for each gene, read the  sp2gene2seq dict and search for each key (the specie) if it contain the gene. If yes, extract the specie id and sequence
    pour chaque gene, parcourir pour chaque clé du meta dico (= pour chaque esp) la cle gene et prendre sp et seq pour ecrire specie2seq
    specie2seq = dict() 
    espwithgene=list()
    for sp in especes_list:
        if gene in sp2gene2seq[sp]:
            specie2seq[sp] = sp2gene2seq[sp][gene]
            espwithgene.append(sp)
    if os.path.exists(PathData + "Genes/{}_no_filter/".format(gene))==False :
        os.mkdir(PathData + "Genes/{}_no_filter/".format(gene))
        sortie = open(PathData + "Genes/{}_no_filter/{}.fasta".format(gene, gene), "w")
        for key, value in specie2seq.items():
            sortie.write(">" + key + "\n" + value + "\n")
        summary = open(PathData + "Genes/summary{}".format(str(split)),"a")
        summary.write("[" + gene + "]" + "[" + str(len(espwithgene)) + "]"  + str(espwithgene) + "\n")

#write a file per species containing the exon position
#its a work by split, have to then manually merge the files

for sp in especes_list:
    if os.path.exists(PathData + "Refexons/{}_no_filter/".format(sp))==False :
        os.mkdir(PathData + "Refexons/{}_no_filter/".format(sp))
    for gene in gene_list:
        gene2ref = dict()
        if gene in sp2gene2ref[sp]:
            gene2ref[gene] = sp2gene2ref[sp][gene]
        for key, value in gene2ref.items():
            sortieref = open(PathData + "Refexons/{}_no_filter/{}.refexons{}_scafinfo".format(sp, sp, str(split)), "a")
            sortieref.write(key + "\t" + sp2gene2scafold[sp][key] + "\t" + sp2gene2scafoldpos[sp][key]+ "\n" + value + "\n")
            #print (key + "\t" + sp2gene2scafold[sp][key] + "\t" + sp2gene2scafoldpos[sp][key]+ "\n" + value + "\n")

