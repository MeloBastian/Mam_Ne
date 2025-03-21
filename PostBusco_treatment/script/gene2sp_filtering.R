a=read.delim("/home/mbastian/data/Enard_postbusco/Genes/genes2sp_summary")
rownames(a) = a$Species_name
a = a[,-1]

library(stringr)
colnames(a) = str_replace_all(colnames(a),"X","")

#rowSums(a)
#colSums(a)

hist(rowSums(a),breaks=50, main="nb genes/sp", xlab="nb genes") #nbr genes per sp
hist(colSums(a),breaks=50, main="nb sp/genes", xlab="nb sp") #nbr sp per genes


#table(colSums(a) > 100)
#table(rowSums(a) > 5000)

#table(table(  sapply(rownames(a),  function(x)  str_split(x,"_")[[1]][1] )  )== 1)


prop_genes_per_species = .7
ncol(a) * prop_genes_per_species
da = a[rowSums(a)/ncol(a) >=  prop_genes_per_species  ,  ] #rm sp without 70% of genes


prop_species_per_genes = .8
prop_species_per_genes * nrow(da)
dg = da[, colSums(da)/nrow(da) >= prop_species_per_genes] #rm genes without 80% sp
hist(rowSums(dg),breaks=50, main="nombre de genes/esp", xlab="nombre genes")
hist(colSums(dg),breaks=50, main="nombre d'esp/genes", xlab="nombre esp")


#write.table(dg, "/home/mbastian/data/Zoonomia_postbusco/filtred_genes_sp", sep="t" )
write.table(dg[1], "/home/mbastian/data/Enard_postbusco/splist_filtred_183sp", sep="t" ) #need manual work
write.table(dg[1,], "/home/mbastian/data/Enard_postbusco/Genes/genelist_filtred_8060genes", sep="\t" )#need manual work

#write.table(rowSums(dg),"/home/mbastian/home/Scripts/nbrgenesbysp")

removesp=a[rowSums(a)/ncol(a)<= prop_genes_per_species , ] #sp removed
write.table(removesp[1], "/home/mbastian/data/Enard_postbusco/spremoved", sep="t" )




