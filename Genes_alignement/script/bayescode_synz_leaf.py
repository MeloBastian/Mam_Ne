gene2count = {}
count_dev = 0
count_all=0
with open("/home/mbastian/data/Enard_postbusco/phylter_vs_bayescode/fulldataset/145sp/genebranchdnds.postmeandev.tab", "r") as file, \
     open("/home/mbastian/data/Enard_postbusco/phylter_vs_bayescode/fulldataset/145sp/bayescode_outlier_synZsup5", "w") as out:

    header = file.readline().split()
    for line in file:
        try:
            fields = line.strip().split()
            if fields[1] == fields[2]:  # only terminal branches
                gene = fields[0]
                count_all+=1
                synz = float(fields[8])
                if synz > 5:  # if deviation in ds is greater than 5
                    count_dev += 1
                    gene2count[gene] = gene2count.get(gene, 0) + 1
                    out.write(gene + "\t" + fields[2] + "\n")
        except IndexError:
            print(f"Skipping line {line!r} due to IndexError")

print(f"Number of sequences with dev > 5: {count_dev}")
print(f"Total number of sequences: {count_all}")


list=gene2count.values()
a=sorted(list)
count=0
for elmt in a:
    if int(elmt)>=5:
        count+=1

