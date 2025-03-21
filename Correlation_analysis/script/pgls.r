library(ape)
#thv=tableTHV_146sp_forR <- read.csv("~/data/Enard_postbusco/THV/tableTHV_146sp_dsom", sep="")

thv=read.table("/home/mbastian/data/Enard_postbusco/THV/6002genes/tableTHV_90sp_6002genes_lowpsallowed_dsom_forR", header=T)
thv=read.table("/home/mbastian/data/Enard_postbusco/THV/6002genes/tableTHV_90sp_6002genes_dsom_forR", header=T)
thv=read.table("/home/mbastian/data/Enard_postbusco/THV/6002genes/tableTHV_144sp_6002genes_lowpsallowed_dsom_forR", header=T)
thv=read.table("/home/mbastian/data/Enard_postbusco/THV/6002genes/tableTHV_144sp_6002genes_dsom_forR", header=T)

phylo =read.tree("/home/mbastian/data/Enard_postbusco/fastcoevol/sept23/1007forphylo.conc.treefile_rooted")
phylo=read.tree("/home/mbastian/data/Enard_postbusco/cuttree10ma/1500genes_subset4.conc.treefile_rooted")
phylo$node.label = NULL

library("caper")
library(nlme)


#regression lineaire
X="maturity"
Y="mass"
subthv <- na.omit(thv[, c("sp",X,Y)])
names(subthv) <- c("sp", "X", "Y")
#thv$pN.pS <- as.numeric(as.character(thv$pN.pS))

modelreg=lm(log(X)~log(Y) , data=subthv)
summary(modelreg)

plot(log(subthv$X)~log(subthv$Y))
abline(modelreg)


#contrastes indep
IC= comparative.data(phylo, thv, sp, na.omit = F , vcv=TRUE )
modelIC=crunch(log(mass)~log(dnds), data=IC)
summary(modelIC) #pgls cor

contrast=caic.table(modelIC) #position des points corrigée apres pgls


plot(contrast$`log(pN.pS`~contrast$`log(dnds)`, xlab="dN/dS", ylab="pN/pS")
plot(contrast$`log(pN.pS`~contrast$`log(dnds)`, xlab="dN/dS", ylab="pN/pS", main="pgls - r=0.21, p-value=0.015")
abline(lm(contrast$`log(pN.pS)`~contrast$`log(dnds)`))

#summary(lm(contrast$`log(dnds)`~contrast$`log(pN.pS)`))

#abline(modelIC)
#summary(contrast$`log(longevity)`)




