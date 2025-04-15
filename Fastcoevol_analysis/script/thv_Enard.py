from statistics import mean

path_THV ="/home/mbastian/data/database"
anagefile =open(path_THV + "/anage.matmasslong", "r")
listspfile=open("/home/mbastian/data/Enard_postbusco/1genus_144splist_E","r")
heterofile=open("/home/mbastian/data/Enard_postbusco/pSwithmask/Enard_pS_masked_GQ150QUAL125_fq02to08_6002genes","r")

specie2mass = dict()
specie2long = dict()
specie2mat = dict()
specie2hetero = dict()
specie2gt = dict()
specie2PS=dict()
specie2pnps=dict()
specie2lentot=dict()
specie2Sobs=dict()

# create a dict per trait with species as keys and "-1" as default value, (correspond to "NA" for fastcoevol) 
listsp = [line.rstrip("\n") for line in listspfile]
for specie in listsp:
    specie2mass[specie] = -1
    specie2long[specie] = -1
    specie2mat[specie] = -1
    specie2hetero[specie] = -1
    specie2gt[specie] = -1
    specie2PS[specie]=-1
    specie2pnps[specie]=-1
    specie2lentot[specie]=-1
    specie2Sobs[specie]=-1



# anage database treatment, create an other dict with keys as species from the db and values correspond to the trait value from the db. 
anagefile.readline()
anage = [line.rstrip("\n") for line in anagefile]

sp2thv = dict()
genre2mass = dict()
genre2mat = dict()
genre2long = dict()
listspnotanage=list()

for elmt in anage:
    listanage = elmt.split(" ")  
    sp2thv[listanage[0]] = [float(listanage[1]), float(listanage[2]), float(listanage[3])]

for sp in listsp:  # write dict with genus as a key and mean traits per genus for the values (if their is multiple sp per genus)
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
        
for sp in listsp:     #want the mean per genus for all species in the analysis
    genremam = sp.split("_")
    if genremam[0] in genre2mass:
        specie2mass[sp] = genre2mass[genremam[0]]
    if genremam[0] in genre2mat:
        specie2mat[sp] = genre2mat[genremam[0]]
    if genremam[0] in genre2long:
        specie2long[sp] = genre2long[genremam[0]]
#print(listspnotanage)


#polymorphism data
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

#generation time computation (not used anymore for the next)
for sp in listsp:
    if specie2mat[sp]!=-1 or specie2mat[sp]!=-1.0:
        if specie2long[sp]!=-1 or specie2long[sp]!=-1.0:
            gt= float(specie2long[sp])*365.25*0.29 + float(specie2mat[sp])
            specie2gt[sp]=gt

# write a traits file formated for fastcoevol
sortie = open("/home/mbastian/data/Enard_postbusco/fastcoevol/dec23/nogt_lowpsallowed/tableTHV_144mam_masked_6002genes_nogt_lowpsallowed", "w")
#sortie = open("/home/mbastian/data/Zoonomia_postbusco/Coevol/juin2022/indeppolym/tableTHV_180Zoonomia_6traits_3648genestail", "w")
sortie.write("#TRAITS" + "\n")
sortie.write(str(len(listsp)) + " 5 mass maturity longevity pS pN/pS" + "\n")

for esp in listsp:
    sortie.write(esp + " " + str(specie2mass[esp]) + " " + str(specie2mat[esp]) + " "  + str(specie2long[esp]) +  " " + str(specie2PS[esp]) + " " + str(specie2pnps[esp]) + "\n")
    print(esp + " " + str(specie2mass[esp]) + " " + str(specie2mat[esp]) + " "  + str(specie2long[esp]) +  " " + str(specie2PS[esp]) + " " + str(specie2pnps[esp]) + "\n")
    #print(esp)
print("done")
