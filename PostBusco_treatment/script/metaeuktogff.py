import sys
name=sys.argv[1]
#metaeuk_file = open("/home/mbastian/data/Enard_postbusco/Refexons/Pongo_abelii_no_filterfullrefexon_scafinfo", "r")
metaeuk_file = open("/beegfs/data/mbastian/Enard_postbusco/Refexons/"+name, "r")




namelite=name.split("_no_filter")
namelite=namelite[0]

output=open("/beegfs/data/mbastian/Enard_postbusco/Refexons/"+namelite+"_phase1.gff", "a")
#output=open("/home/mbastian/data/Zoonomia_postbusco/Refexons/"+namelite+"crochet.gff", "w")
#metaeuk_file = open("/home/mbastian/data/Zoonomia_postbusco/Refexons/Refexons_V2/"+name, "r")
metaeuk=metaeuk_file.readlines()
i=0
linetoread=0
while i<len(metaeuk):
    info=metaeuk[linetoread]
    linetoread+=1
    i+=1
    info=info.split()
    #print(info)
    scafold = info[1]
    gene=info[0]
    ref=metaeuk[linetoread]
    linetoread += 2 #bc \n
    i += 2
    ref=ref.split("|")
    exons=ref[3:]
    phase=0
    for elmt in exons: #(10[10]:30[30]:20[20])
        elmt=elmt.split(":") # (10[10]),(30[30]),(20[20])
        #lower=elmt[0].split("[")[0] #outbracket
        lower = elmt[0].split("[")[1][:-1] #inbracket
        #high=elmt[1].split("[")[0]
        high = elmt[1].split("[")[1][:-1]
        strand ="N"
        start="N"
        end="N"

        if int(lower)<int(high):
            strand="+"
            start=(int(lower)+1)
            end=(int(high)+1)

        else:
            strand = "-"
            start=(int(high)+1)
            end=(int(lower)+1)

        output.write(scafold+ "\t"+ "CDS" + "\t"+ str(start) + "\t" + str(end) +"\t"+ strand + "\t" + str(phase) + "\t" + gene + "\n")
        #print(scafold+ "\t"+ "CDS" + "\t"+ start + "\t" + end +"\t"+ strand + "\t" + str(phase) + "\t" + gene + "\n")

        rest = (int(end) - int(start) + 1 - phase) % 3
        if rest == 2 :
            phase=1
            print(scafold+ "\t"+ str(start) + "\t" + str(end) +"\t" + gene )
            print("changement en phase 1")
        if rest == 1 :
            phase =2
            print(scafold+ "\t"+ str(start) + "\t" + str(end) +"\t" + gene )
            print("changement en phase 2")
        if rest == 0:
            phase =0
print("done " + namelite)
