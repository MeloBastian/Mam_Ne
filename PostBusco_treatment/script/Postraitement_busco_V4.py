
import os
from Bio import SeqIO
from Bio.Seq import Seq
from ali import translate_nucseq
import sys

split = sys.argv[1]
#split = 0
print("split: " + split)
PathBanque = sys.argv[2]
PathData = sys.argv[3]


#especes_list_path = PathBanque + "Zoonomia_Busco/Busco_splited_species/listmam"
especes_list_path = sys.argv[4]
text_file = open(especes_list_path, "r")
lines = text_file.readlines()
lines = [f.replace(".fna\n", "") for f in lines]
especes_list = [f.replace("mam_", "") for f in lines]

gene_list = list()



#especes_list=especes_list[238:]


print(especes_list)

#especes_list=("Hystrix_cristata")



sp2gene2seq = dict()  # chaque cle espece contient un dico gene2seq
sp2gene2ref = dict()  # chaque cle esp contient un dico gene2ref
sp2gene2scafold = dict()
sp2gene2scafoldpos = dict()

##########################################################################################################################
for sp in especes_list:  # pour chaque sp, faire un gene2seq de toutes les sequences nuc apres filtrage et un gene2ref pour recuperer les infos exons
    print()
    print(sp)
    genome_dict = SeqIO.to_dict(SeqIO.parse("/beegfs/banque/mbastian/Enard/genomes/"+"mam_" + sp + ".fna", "fasta"))
    path_sp = PathBanque + "/Busco_splited_species/mam_" + sp + ".fna/"
    gene2seq = dict()  # va contenir toutes les sequences nuc apres filtrage
    gene2ref = dict()  # va recuperer les info exons pour chaque gene
    gene2count = dict()  # compte l'occurence de chaque gene suppose unique dans les differents splits d'une meme esp
    gene2aaseq = dict()  # va contenir toutes les seq aa pour une sp
    listegenesinglesp = list()  # va contenir les identifiants de genes qui sont reellements single copie dans tout le genome
    gene2scafold = dict() #for the scafold ref
    gene2scafoldpos = dict() #for the position in the scafold

    compte_seqdiffrerun = 0
    compte_seqdiffinitial = 0
    compte_taille = 0
    compte_falsesingle = 0
    compte_exon=0
    compte_exons_out_false=0
    compte_exons_in_false=0
    nocorresp=0
    nbrseq=0
    nbrwrongseq_in=0
    nbrwrongseq_out=0
    # premier parcours des splits de l'especes afin d'etablir listegenesinglesp et gene2aaseq
    #singleseqpath="/home/mbastian/Bureau/busco_mam_Hystrix_cristata.fna_0/run_mammalia0_odb10/busco_sequences/single_copy_busco_sequences/"
    singleseqpath = path_sp + "busco_mam_" + sp + ".fna_" + str(split) + "/run_mammalia" + str(split) + "_odb10/busco_sequences/single_copy_busco_sequences/"  # contient un fichier par gene annonce uniseq
    for gene in os.listdir(singleseqpath):
        id = gene[:-4]  # supprime le ".faa"
        if id not in gene2count:  # prend uniquement la premiere occurence de la seq aa (normalement toutes id)
            gene2count[id] = 1
            aapath = singleseqpath + id + ".faa"  # permet de n'ouvrir que les fichiers de cette boucle
            if os.path.isfile(aapath):  # parfois le file est exotique...
                file = open(aapath, "r")
                record = SeqIO.read(file, "fasta")
                idseq =str(record.id)
                idseq =idseq.split(":")
                idscafold=idseq[0]
                idscafold_pos=idseq[1]
                #print (idscafold + idscafold_pos)

                gene2scafold[id]=idscafold
                gene2scafoldpos[id]=idscafold_pos
                gene2aaseq[id] = str(record.seq)
        else:
            gene2count[id] += 1
            if gene2count[id] >1 :
                print("gene" +id+" non unique")
    for keys in gene2count.keys():  # on recupere ceux qui n'ont ete observes qu'une fois: valeur de comptage =1
        if gene2count[keys] == 1:
            listegenesinglesp.append(keys)

    # second parcours des splits de l'espece pour recuperer les sequences codon des genes de listegenesinglesp et etablir plusiuers controles qualite

    singleseqpath = path_sp + str(split) + "/run_mammalia_odb10/busco_sequences/single_copy_busco_sequences/"
    fichier = "mam_" + sp + ".fna" + ".codon.fas"  # lecture du fichier initial/rerun codon pour chaque split, doit enlever denomination "busco"
    #pathinitial="/home/mbastian/Bureau/busco_mam_Hystrix_cristata.fna_0/run_mammalia0_odb10/metaeuk_output/initial_results/mam_Hystrix_cristata.fna.codon.fas"
    pathinitial = path_sp + "busco_mam_" + sp + ".fna_" + str(split) + "/run_mammalia" + str(split) + "_odb10/metaeuk_output/initial_results/" + fichier
    #pathrerun="/home/mbastian/Bureau/busco_mam_Hystrix_cristata.fna_0/run_mammalia0_odb10/metaeuk_output/rerun_results/mam_Hystrix_cristata.fna.codon.fas"
    pathrerun = path_sp + "busco_mam_" + sp + ".fna_" + str(split) + "/run_mammalia" + str(split) + "_odb10/metaeuk_output/rerun_results/" + fichier

    # traitement des sorties initialresults
    alreadytreated=list()
    if os.path.exists(pathinitial):  # verifie si le fichier initial codon exist
        initialcodon = open(pathinitial, "r")
        while 1:
            idfull = initialcodon.readline()  # lis la premiere ligne : info exon
            if idfull == "":  # quand aura atteint la derniere ligne
                break
            seq = initialcodon.readline().strip()  # lis la deuxieme ligne : sequence
            idlist = idfull.split("_")
            id = idlist[0][1:]  # id du gene sans chevron
            id2 = id.split("|")  # ?
            id3 = id2[0]  # ?
            reffull = idfull.split("|", maxsplit=5)
            scafinitial=reffull[1] #sc  af ref from initial files
            ref = reffull[-1]  # prend uniquement info exon
            #print(ref)


            if id3 in listegenesinglesp and id3 not in alreadytreated :  # si le gene est reellement single copy et que la bonne seq n'as pas deja ete trouvee
                nbrseq += 1

                # the CDS have to be multiple of 3
                if len(seq) % 3 != 0:
                    print("la sequence n'est pas un multiple de 3", sp, id3)
                    compte_taille += 1
                else:
                    seq_up = seq.upper()
                    # the trad of the sequence codon in initial results have to be the same than the seq aa
                    translate_nuc = translate_nucseq(seq_up)
                    translate_nuc_up = translate_nuc.upper()
                    aaseq = gene2aaseq[id3]
                    aaseq_up = aaseq.upper()
                    if translate_nuc_up != aaseq_up:
                        # print("sequence nuc", sp, id3, "differentes de la sequence aa")
                        compte_seqdiffinitial += 1

                    # work on the exons from CDS and from genome
                    scafold = gene2scafold[id3]
                    if scafold == scafinitial : #test if the scafold id from busco seq are the same than in the initial fiels
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
                            end_in = int(end_in)   # +1 to have the last position incle in the interval in python

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

                        if errorinfullseq_in != 0:  # after the rewrite of the CDS, if one exon is not well
                            nbrwrongseq_in += 1
                        else:
                            if CDSfromgenome != seq_up:  # try if the CDS reconstituated from the genome is the same than the one from busco ("!=" for now but "not in" before)
                                print("CDS from genome different from Busco")
                                # nbrwrongCDS+=1
                            else:
                                gene2seq[id3] = seq  # si tous les filtres sont passes, recuperer la sequence dans gene2seq
                                gene2ref[id3] = ref
                                alreadytreated.append(id3)  # if a good CDS is found with all the filter corect, do not try to find an other CDS for this gene
                    # else:
                    #     print("id scafold from initial files is different than busco sequence file ")

    else:
        print(pathinitial+ "doesn't exist")  # si le fichier initial codon n'existe pas
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

    # traitement des sorties rerunresults
    if os.path.exists(pathrerun):  # verifie si le fichier initial codon exist
        initialcodon = open(pathrerun, "r")
        while 1:
            idfull = initialcodon.readline()  # lis la premiere ligne : info exon
            if idfull == "":  # quand aura atteint la derniere ligne
                break
            seq = initialcodon.readline().strip()  # lis la deuxieme ligne : sequence
            idlist = idfull.split("_")
            id = idlist[0][1:]  # id du gene sans chevron
            id2 = id.split("|")  # ?
            id3 = id2[0]  # ?
            reffull = idfull.split("|", maxsplit=5)
            ref = reffull[-1]  # prend uniquement info exon
            #print(ref)


            if id3 in listegenesinglesp and id3 not in alreadytreated :  # si le gene est reellement single copy
                nbrseq += 1

                # the CDS have to be multiple of 3
                if len(seq) % 3 != 0:
                    print("la sequence n'est pas un multiple de 3", sp, id3)
                    compte_taille += 1
                else:
                    seq_up = seq.upper()
                    # the trad of the sequence codon in initial results have to be the same than the seq aa
                    translate_nuc = translate_nucseq(seq_up)
                    translate_nuc_up = translate_nuc.upper()
                    aaseq = gene2aaseq[id3]
                    aaseq_up = aaseq.upper()
                    if translate_nuc_up != aaseq_up:
                        # print("sequence nuc", sp, id3, "differentes de la sequence aa")
                        compte_seqdiffinitial += 1

                    # work on the exons from CDS and from genome
                    scafold = gene2scafold[id3]
                    if scafold == scafinitial : #test if the scafold id from busco seq are the same than in the initial fiels
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
                            end_in = int(end_in)

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

                            if len(seq_exon)!=0: #some exons look like having a size 0
                                if str(seq_exon).upper() not in str(seq_up):
                                    compte_exons_in_false += 1
                                    errorinfullseq_in += 1
                                    nocorresp+=1
                                    print("genome seq_exon not in busco CDS ")

                                else:
                                    CDSfromgenome += str(seq_exon)  # concatenate the exons
                            else:
                                print("len exon = 0") #be careful to the previous +1
                                compte_exons_in_false += 1
                                errorinfullseq_in += 1

                        if errorinfullseq_in != 0:  # after the rewrite of the CDS, if one exon is not well
                            nbrwrongseq_in += 1
                        else:
                            if CDSfromgenome != seq_up:  # try if the CDS reconstituated from the genome is the same than the one from busco ("!=" for now but "not in" before)
                                print("CDS from genome different from Busco")
                                # nbrwrongCDS+=1
                            else:
                                gene2seq[id3] = seq  # si tous les filtres sont passes, recuperer la sequence dans gene2seq
                                gene2ref[id3] = ref
                                alreadytreated.append(id3)  # if a good CDS is found with all the filter corect, do not try to find an other CDS for this gene
                    # else:
                    #     print("id scafold from rerun files is different than busco sequence file ")

    else:
        print(pathinitial+ "doesn't exist")  # si le fichier initial codon n'existe pas
    nbrexonfalsein=(100*compte_exons_in_false)/compte_exon
    tauxwronseq_in=(100*nbrwrongseq_in)/nbrseq
    print("initial+rerun result : % exon false "+str(nbrexonfalsein))
    print("initial+rerunn result : % CDS false  "+str(tauxwronseq_in))
    print("no coresp exon from genome and CDS from busco", str(nocorresp))



    print(sp, "nombre de genes single pré filtrage =", len(listegenesinglesp),", nombre de sequences gardées =", len(gene2seq), " --> taux de sequences gardées =",(len(gene2seq) / len(listegenesinglesp)) * 100, "%")






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
# faire une liste de gene
for gene in gene_list:  # pour chaque gene, parcourir pour chaque clé du meta dico (= pour chaque esp) la cle gene et prendre sp et seq pour ecrire specie2seq
    specie2seq = dict()  # va contenir pour un gene le dico esp:seq
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

        # probleme car travail par split alors que je veux un unique fichier par esp pour les refexons --> a concatener plus tard
