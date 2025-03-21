from statistics import mean

path_THV ="/home/mbastian/data/database"
anagefile =open(path_THV + "/anage.matmasslong", "r")
#genomesize= open(path_THV + "/genomesizelite.txt","r")
listspfile=open("/home/mbastian/data/Enard_postbusco/1genus_144splist_E","r")
heterofile=open("/home/mbastian/data/Enard_postbusco/pSwithmask/Enard_pS_masked_GQ150QUAL125_fq02to08_6002genes","r")
#heterofile=open("/home/mbastian/data/Zoonomia_postbusco/Coevol/juin2022/indeppolym/tablepnps_fullinfo_notali_3648genestail","r")
#listspfile=open("/home/mbastian/data//Orthomam/listsp_Orthomam_lite","r")
#heterofile=open("/home/bastian/Bureau/zoonomia/table_hetero.txt", "r")

specie2mass = dict()
specie2long = dict()
specie2mat = dict()
specie2hetero = dict()
specie2gt = dict()
#specie2GS = dict()
specie2PS=dict()
specie2pnps=dict()
specie2lentot=dict()
specie2Sobs=dict()

# remplir les cles des dico
listsp = [line.rstrip("\n") for line in listspfile]
#listsp=[line[4:-4] for line in listsp]

for specie in listsp:
    specie2mass[specie] = -1
    specie2long[specie] = -1
    specie2mat[specie] = -1
    specie2hetero[specie] = -1
    specie2gt[specie] = -1
    #specie2GS[specie] = -1
    specie2PS[specie]=-1
    specie2pnps[specie]=-1
    specie2lentot[specie]=-1
    specie2Sobs[specie]=-1



# traitement bd anage, a du supprimer les tabulations dans le fichier
anagefile.readline()
anage = [line.rstrip("\n") for line in anagefile]

sp2thv = dict()
genre2mass = dict()
genre2mat = dict()
genre2long = dict()
listspnotanage=list()
for elmt in anage:
    listanage = elmt.split(" ")  # nom,  maturite, masse, longevite
    sp2thv[listanage[0]] = [float(listanage[1]), float(listanage[2]), float(listanage[3])]


for sp in listsp:  # creation dico avec moyenne par genre pour les 3 thv
    genremam = sp.split("_") #(genre, sp)
    genremass = list()
    genremat = list()
    genrelong = list()
    for key, thv in sp2thv.items():
        genreanage = key.split("_")
        if genremam[0] == genreanage[0]:
            if thv[0] != -1:
                genremat.append(thv[0])
            if thv[1] != -1:
                genremass.append(thv[1])
            if thv[2] != -1:
                genrelong.append(thv[2])
    if genremass:
        genre2mass[genremam[0]] = mean(genremass)
    if genremat:
        genre2mat[genremam[0]] = mean(genremat)
    if genrelong:
        genre2long[genremam[0]] = mean(genrelong)
        
for sp in listsp:     #want the mean per genus for all
    # if sp in sp2thv:
    #     specie2mat[sp]=sp2thv[sp][0]
    #     specie2mass[sp] = sp2thv[sp][1]
    #     specie2long[sp] = sp2thv[sp][2]
    # else:
    #     listspnotanage.append(sp)
    genremam = sp.split("_")
    if genremam[0] in genre2mass:
        specie2mass[sp] = genre2mass[genremam[0]]
    if genremam[0] in genre2mat:
        specie2mat[sp] = genre2mat[genremam[0]]
    if genremam[0] in genre2long:
        specie2long[sp] = genre2long[genremam[0]]
#print(listspnotanage)


#traitement hetero
#sptoremove=["Sigmodon_hispidus","Muscardinus_avellanarius", "Acomys_cahirinus", "Cheirogaleus_medius","Przewalskium_albirostris","Litocranius_walleri","Alouatta_palliata", "Mastomys_coucha", "Sousa_chinensis","Diceros_bicornis", "Beatragus_hunteri", "Cephalophus_harveyi"]
sptoremove=["Acomys_cahirinus", "Cheirogaleus_medius","Przewalskium_albirostris","Litocranius_walleri", "Mastomys_coucha",  "Cephalophus_harveyi"]

hetero = [line.rstrip("\n") for line in heterofile]
for elmt in hetero[1:]:
   listhetero= elmt.split("\t")
   if listhetero[0]+"_E" in listsp and listhetero[0] not in sptoremove:
       specie2PS[listhetero[0]+"_E"]=listhetero[1]
       specie2pnps[listhetero[0]+"_E"]=listhetero[3]
       specie2lentot[listhetero[0]+"_E"]=listhetero[6]
       #specie2Sobs[listhetero[0]]=listhetero[7]

#temps generation
for sp in listsp:
    if specie2mat[sp]!=-1 or specie2mat[sp]!=-1.0:
        if specie2long[sp]!=-1 or specie2long[sp]!=-1.0:
            gt= float(specie2long[sp])*365.25*0.29 + float(specie2mat[sp])
            specie2gt[sp]=gt

#genomesize
# listGSprst=list()
# GS = [line.rstrip("\n") for line in genomesize]
# for sp in listsp:
#     for elmt in GS:
#         listGS=elmt.split("\t")
#         if listGS[0] == sp:
#             specie2GS[listGS[0]]=float(listGS[1])
#             listGSprst.append(sp)
# listGSprst=list(set(listGSprst))
# print("listgsprst : ", listGSprst)
# print("nombre de GS:", len(listGSprst))

# ecrire dans un fichier pour coevol
sortie = open("/home/mbastian/data/Enard_postbusco/fastcoevol/dec23/nogt_lowpsallowed/tableTHV_144mam_masked_6002genes_nogt_lowpsallowed", "w")
#sortie = open("/home/mbastian/data/Zoonomia_postbusco/Coevol/juin2022/indeppolym/tableTHV_180Zoonomia_6traits_3648genestail", "w")
sortie.write("#TRAITS" + "\n")
sortie.write(str(len(listsp)) + " 5 mass maturity longevity pS pN/pS" + "\n")

for esp in listsp:
    #sortie.write("mam_" + esp + ".fna" + " " + str(specie2mass[esp]) + " " + str(specie2mat[esp]) + " "  + str(specie2long[esp]) + " " +  str(specie2gt[esp]) + " " + str(specie2GS[esp]) + " " + str(specie2hetero[esp]) + "\n")
    sortie.write(esp + " " + str(specie2mass[esp]) + " " + str(specie2mat[esp]) + " "  + str(specie2long[esp]) +  " " + str(specie2PS[esp]) + " " + str(specie2pnps[esp]) + "\n")
    print(esp + " " + str(specie2mass[esp]) + " " + str(specie2mat[esp]) + " "  + str(specie2long[esp]) +  " " + str(specie2PS[esp]) + " " + str(specie2pnps[esp]) + "\n")
    #print(esp)
print("done")
