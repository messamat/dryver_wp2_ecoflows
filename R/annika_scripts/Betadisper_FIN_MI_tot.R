library(vegan)
inv.spec= read.csv2("D:/Users/mykrah/Omat tiedostot/Projektit/DRYVER/WP 2/Metacommunity modeling/R files/FIN_MI_6.csv", header=T,row.names=1, sep=";", dec=".")

# maakohtaisesti eli tee loop ? 

groups<-factor(c(rep(1,20), rep(2,21), rep(3,21), rep(4, 18), rep(5,21), rep(6,21)), labels=c("Apr", "May", "Jun", "Aug", "Sep", "Dec"))

#Function for Jaccard
jac =function(x) vegdist(x, "jaccard", binary =TRUE)
#Null model
nul <- oecosimu(resp, jac, nsimul = 999, "quasiswap")

##STANDARDIZED EFFECT SIZE

dis.jac.z <- vegdist(resp)
dis.jac.z[] <- nul$oecosimu$z

disdif.jac.z<-as.dist((dis.jac.z-min(dis.jac.z))/(max(dis.jac.z)-min(dis.jac.z))) ### standar



###BETADISPER#####Vertaa kampanjoiden välisiä eroja drn:n sisällä
beta.dis.jac.z <-betadisper(disdif.jac.z, groups) # exp$Campaign
beta.dis.jac.z
anova(beta.dis.jac.z)
permutest(beta.dis.jac.z, pairwise = TRUE)
plot(beta.dis.jac.z)
boxplot(beta.dis.jac.z,cex.lab=1.1 )
TukeyHSD(beta.dis.jac.z)
## Next step is saving the distances results --> remember to name the file e.g. here: distances.csv
write.csv2(as.matrix(beta.dis.jac.z$distances), file="M:/gDRYvER/WP 2/Metacommunity models/Finland/distances_tot2.csv")


######################### ilman nollamallia (havaitulla jaccardin indeksillä eli ei ole tehty mitään datalle)
jac.obs <-vegdist(inv.spec, "jaccard", binary =TRUE)
beta.jac.obs <-betadisper(jac.obs, groups)
beta.jac.obs
anova(beta.jac.obs)
permutest(beta.jac.obs, pairwise = TRUE)
plot(beta.jac.obs)
boxplot(beta.jac.obs)









################WITHIN GROUP NULL MODEL############# tätä voisi käyttää paikkakohtaisen (n = 6) dissimilariteetin laskemiseen
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
oecosimu(inv.spec, foo, "quasiswap", nsimul = 999, groups = groups)


fi.taxa= read.csv2("D:/Users/mykrah/Omat tiedostot/Projektit/DRYVER/WP 2/Metacommunity modeling/R files/FIN6_taxa.csv", header=T,row.names=1, sep=";", dec=".")