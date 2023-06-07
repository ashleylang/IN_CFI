#------ AMF 1 -----------

getwd()

setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/AMF")
myfileOTU <- read.csv("Seq_1.csv", head =TRUE, row.names=1)
myfileOTU <- t(myfileOTU)
 
PlantDiv <- read.csv("PD_1.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PS_1.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PP_1.csv", head =TRUE, row.names=1)
Soils <- read.csv("Fe_1.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_1.csv", head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

#determine MEMs

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,1:8)
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)
dumdum<-as.matrix(dummy)


#run full matrix model to determine significant components


dbrda <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils  + finalGEO, distance = "bray")


mardb <- Anova(dbrda, type = 3)

namemar <- subset(mardb, mardb[,3] <= 0.050)

signamesmar <- as.factor(row.names(namemar))

print(signamesmat)

#run reduced matrix model to determine correct p-values

redmodmar<-dbrda(ddat ~ soils + finalGEO)


aovredmodmar <- Anova(redmodmar,type = 3)



write.csv(aovredmodmar, "AMF_1mar.csv")



#run full variable model to determine significant components

dbrda<-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)

mardb <- Anova(dbrda, type = 3)



namemar <- subset(mardb, mardb[,3] <= 0.050)

signamesmar <- as.factor(row.names(namemar))

print(signamesmar)

redmodmar<-dbrda(ddat ~Soils$cn +Soils$pH+finalGEO)
aovredmodmar <- Anova(redmodmar, type = 3)


nameter <- subset(aovredmodmar, aovredmodmar[,4] <= 0.050)

signamester <- as.factor(row.names(nameter))

print(signamester)



#dataframe with all models
AMF1 <- data.frame(maraovredmodmat,  maraovredmodvar)


#csv for easy copypasta
write.csv(AMF1, "AMF_0_5_chi.csv")

VPAMF <- varpart(ddat,  soils, dumdum)
VPAMF

plot(VPAMF)
windows()
bac1varpart <- 
  plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))



#------ AMF 2---------------
rm(list=ls())
# myfileOTU <- read.csv(file.choose(), head =TRUE, row.names=1)
# myfileOTU<-t(myfileOTU)

setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/AMF")


myfileOTU <- read.csv("IndianaCFIAMF_seqtab_prop_5_10_2022_02_24.csv", head =TRUE, row.names=1)
myfileOTU<-t(myfileOTU)


PlantDiv <- read.csv("PD_2.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PS_2.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PP_2.csv", head =TRUE, row.names=1)
Soils <- read.csv("Fe_2.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_2.csv", head =TRUE, row.names=1)
#myfileENV<-read.csv(file.choose(), head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,c(4)])
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)


#run full matrix model to determine significant components


dbrda <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils  + finalGEO, distance = "bray")


mardb <- Anova(dbrda, type = 3)

namemar <- subset(mardb, mardb[,3] <= 0.050)

signamesmar <- as.factor(row.names(namemar))

print(signamesmar)

#run reduced matrix model to determine correct p-values

redmodmar<-dbrda(ddat ~ plantdiv + plantProd + soils + finalGEO)


aovredmodmar <- Anova(redmodmar,type = 3)



write.csv(aovredmodmar, "AMF_1mar.csv")



#run full variable model to determine significant components

dbrda<-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)

mardb <- Anova(dbrda, type = 3)



namemar <- subset(mardb, mardb[,3] <= 0.050)

signamesmar <- as.factor(row.names(namemar))

print(signamesmar)

redmodmar<-dbrda(ddat ~PlantStructure$VAI  + PlantDiv$AMdominance + PlantProd$stand_age_years + Soils$cn +Soils$pH+  Soils$Fe_percent +finalGEO)
aovredmodmar <- Anova(redmodmar, type = 3)


nameter <- subset(aovredmodmar, aovredmodmar[,3] <= 0.050)

signamester <- as.factor(row.names(nameter))

print(signamester)



#dataframe with all models
AMF2 <- data.frame(maraovredmodmat,  maraovredmodvar)


#csv for easy copypasta
write.csv(AMF2, "AMF_5_10_chi.csv")

VPAMF2 <- varpart(ddat,  PlantDiv, PlantProd, soils)
VPAMF2

plot(VPAMF)
windows()
bac1varpart <- plot (VPAMF2, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))



#-------------Bacteria 1---------------

rm(list=ls())
 
setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/Bacteria")
myfileOTU <- read.csv("IndianaCFI16S_Prop_seqtab_0_5_2021_04_16.csv", head =TRUE, row.names=1)
myfileOTU<-t(myfileOTU)
# 
# setwd("C:/Users/Earth-Person/Desktop/Lang_etal_data_archive/data/betadiversity/AMF")
# 
# PlantDiv <- read.csv(file.choose(), head =TRUE, row.names=1)
# PlantStructure <- read.csv(file.choose(), head =TRUE, row.names=1)
# PlantProd <- read.csv(file.choose(), head =TRUE, row.names=1)
# Soils <- read.csv(file.choose(), head =TRUE, row.names=1)
# Geo <- read.csv(file.choose(), head =TRUE, row.names=1)
# #myfileENV<-read.csv(file.choose(), head =TRUE, row.names=1)

PlantDiv <- read.csv("PlantTreeDiv_0_5_02_2022.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PlantStructural_0_5_02_2022.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PlantProd_0_5_02_2022.csv", head =TRUE, row.names=1)
Soils <- read.csv("fe_pH_0_5_2022_02_15.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_0_5_02_2022.csv", head =TRUE, row.names=1)
#myfileENV<-read.csv(file.choose(), head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,c(4)])
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)




dbrdafullmat <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils + finalGEO)


mardbrdafullmat <- Anova(dbrdafullmat, type = 3)


namemar <- as.factor(row.names(subset(mardbrdafullmat, mardbrdafullmat[,3] <= 0.10)))

print(namemar)


#reduced marginal model based on sig value

marredmodmat<-dbrda(ddat ~ plantStructure + plantdiv    +   soils    + finalGEO )

maraovredmodmat <- Anova(marredmodmat,type = 3)



# #models for individual variables

dbrdafullvar <-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)


marfullvar <- Anova(dbrdafullvar, type = 3)

#names for marginal model

namemar <- as.factor(row.names(subset(marfullvar, marfullvar[,3] <= 0.10)))

print(namemar)


#reduced marginal model based on sig value

marredmodvar<-dbrda(ddat ~ PlantDiv$AMdominance + Soils$cn  +  Soils$pH   + finalGEO)
maraovredmodvar <- Anova(marredmodvar,type = 3)



#this step will tell you which of the five matrices to include in the varpart
VPAMF <- varpart(ddat,  plantStructure,plantdiv,  soils)
VPAMF

plot(VPAMF)
windows()
bac1varpart <- plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))

#-----------Bacteria 2 ----------
rm(list=ls())
setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/Bacteria")
myfileOTU <- read.csv("IndianaCFI16S_Prop_seqtab_5_10_2022_02_16.csv", head =TRUE, row.names=1)
myfileOTU<-t(myfileOTU)

PlantDiv <- read.csv("PlantTreeDiv__5_10_02_2022.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PlantStructural_5_10_02_2022.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PlantProd_5_10_02_2022.csv", head =TRUE, row.names=1)
Soils <- read.csv("fe_pH_5_10_2022_02_15.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_5_10_02_2022.csv", head =TRUE, row.names=1)
#myfileENV<-read.csv(file.choose(), head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,c(4)])
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)




dbrdafullmat <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils + finalGEO)

mardbrdafullmat <- Anova(dbrdafullmat,type = 3)

#some script to print the names
namemar <- as.factor(row.names(subset(mardbrdafullmat, mardbrdafullmat[,3] <= 0.10)))

print(namemar)


#reduced marginal model based on sig value

marredmodmat<-dbrda(ddat ~ plantStructure + plantdiv    +   soils    + finalGEO )

maraovredmodmat <- Anova(marredmodmat,type = 3)



VPAMF <- varpart(ddat,  plantStructure,plantdiv,  soils)
VPAMF

plot(VPAMF)

bac1varpart <- plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))

windows()
bac2varpart <- plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))

# #models for individual variables

dbrdafullvar <-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)


marfullvar <- Anova(dbrdafullvar, type = 3)

#names for marginal model

namemar <- as.factor(row.names(subset(marfullvar, marfullvar[,3] <= 0.10)))

print(namemar)


#reduced marginal model based on sig value

marredmodvar<-dbrda(ddat ~ PlantDiv$AMdominance +  Soils$cn   + Soils$pH   + finalGEO)
maraovredmodvar <- Anova(marredmodvar,type = 3)

 
maraovredmodvar$row.names <- row.names(maraovredmodvar)

maraovredmodmat$row.names <- row.names(maraovredmodmat)

#dataframe with all models
bacteria1 <- data.frame(maraovredmodmat)

bacteria1[0:2,5:9] <- maraovredmodvar






#-------------------- Fungi 1----------------
rm(list=ls())
setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/Fungi")
myfileOTU <- read.csv("IndianaCFIITSProp_seqtab_0_5_2022_02_24_1.csv", head =TRUE, row.names=1)
myfileOTU<-t(myfileOTU)

PlantDiv <- read.csv("PlantTreeDiv_0_5_02_2022.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PlantStructural_0_5_02_2022.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PlantProd_0_5_02_2022.csv", head =TRUE, row.names=1)
Soils <- read.csv("fe_pH_0_5_2022_02_15.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_0_5_02_2022.csv", head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,c(4)])
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)



#run full matrix model to determine significant components

dbrdafullmat <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils + finalGEO)


#anova models

#marginal
mardbrdafullmat <- Anova(dbrdafullmat, type = 3)



#some script to print the names
namemar <- as.factor(row.names(subset(mardbrdafullmat, mardbrdafullmat[,3] <= 0.10)))

print(namemar)



VPAMF <- varpart(ddat,  plantStructure,plantdiv,  soils)
#VPAMF

plot(VPAMF)

plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))


fun1varpart <- plot(VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))

windows()
bac2varpart <- plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'Soils'), bg = c('navy', 'tomato', 'goldenrod'))

#reduced marginal model based on sig value

marredmodmat<-dbrda(ddat ~ plantStructure + plantdiv    +   soils   + finalGEO  )

maraovredmodmat <- Anova(marredmodmat,type = 3)


# #models for individual variables

dbrdafullvar <-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)


marfullvar <- Anova(dbrdafullvar, type = 3)



#names for marginal model

namemar <- as.factor(row.names(subset(marfullvar, marfullvar[,3] <= 0.10)))

print(namemar)






#reduced marginal model based on sig value

marredmodvar<-dbrda(ddat ~ PlantStructure$VCI + PlantDiv$AMdominance + Soils$cn    +         Soils$pH + Soils$Fe_percent      + finalGEO )
maraovredmodvar <- Anova(marredmodvar,type = 3)



maraovredmodvar$row.names <- row.names(maraovredmodvar)

maraovredmodmat$row.names <- row.names(maraovredmodmat)

#dataframe with all models
fungi1 <- data.frame(maraovredmodmat,  maraovredmodvar)



#csv for easy copypasta
write.csv(fungi1, "Fungi_0_5_chi.csv")




#-------------------- Fungi 2----------------
rm(list=ls())
setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/Fungi")
myfileOTU <- read.csv("IndianaCFIITSProp_seqtab_5_10_2022_02_24_1.csv", head =TRUE, row.names=1)
myfileOTU<-t(myfileOTU)

PlantDiv <- read.csv("PlantTreeDiv_5_10_02_2022.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PlantStructural_5_10_02_2022.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PlantProd_5_10_02_2022.csv", head =TRUE, row.names=1)
Soils <- read.csv("fe_pH_5_10_2022_02_15.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_5_10_02_2022.csv", head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,c(4)])
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)

#run full matrix model to determine significant components

dbrdafullmat <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils + finalGEO)


#anova models

#marginal
mardbrdafullmat <- Anova(dbrdafullmat, type = 3)




#some script to print the names
namemar <- as.factor(row.names(subset(mardbrdafullmat, mardbrdafullmat[,3] <= 0.10)))

print(namemar)





VPAMF <- varpart(ddat,  plantdiv,  soils)
#VPAMF

plot(VPAMF)

plot (VPAMF, digits = 3, Xnames = c( 'Plant Diversity', 'Soils'), bg = c( 'tomato', 'goldenrod'))





#reduced marginal model based on sig value

marredmodmat<-dbrda(ddat ~  plantdiv    +   soils   + finalGEO  )

maraovredmodmat <- Anova(marredmodmat,type = 3)

# 
# #models for individual variables
# 
# 

dbrdafullvar <-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)


marfullvar <- Anova(dbrdafullvar, type = 3)



#names for marginal model

namemar <- as.factor(row.names(subset(marfullvar, marfullvar[,3] <= 0.10)))

print(namemar)


#reduced marginal model based on sig value

marredmodvar<-dbrda(ddat ~    PlantDiv$AMdominance +  Soils$cn    +         Soils$pH   + Soils$Fe_percent + finalGEO)
maraovredmodvar <- Anova(marredmodvar,type = 3)


#names for sequential model
# 
# nameter <- as.factor(row.names(subset(terfullvar, terfullvar[,4] <= 0.050)))
# 
# print(nameter)
# 
# #reduced sequential model based on individual variables
# 
# terredmodvar<-dbrda(ddat ~ PlantStructure$vert.sd + PlantStructure$VCI   +  PlantDiv$AMdominance+   Soils$cn       +        Soils$pH         +      Soils$Fe_percent    + Condition(finalGEO))
# aovredmodter <- Anova(redmodter,by="terms", permu=999)
# 

maraovredmodvar$row.names <- row.names(maraovredmodvar)

maraovredmodmat$row.names <- row.names(maraovredmodmat)

#dataframe with all models
fungi2 <- data.frame(maraovredmodmat,  maraovredmodvar)


#csv for easy copypasta
write.csv(fungi2, "Fungi_5_10_chi.csv")





#-------------------- EMF 1----------------

setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/EMF")
myfileOTU <- read.csv("IndianaCFIITSProp_seqtabEcto_0_5_2022_02_24_1.csv", head =TRUE, row.names=1)
myfileOTU<-t(myfileOTU)

PlantDiv <- read.csv("PlantTreeDiv_0_5_02_2022_1.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PlantStructural_0_5_02_2022_1.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PlantProd_0_5_02_2022_1.csv", head =TRUE, row.names=1)
Soils <- read.csv("fe_pH_0_5_2022_02_15_1.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_0_5_02_2022_1.csv", head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,c(4)])
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)

#run full matrix model to determine significant components


dbrdafullmat <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils )


#anova models

#marginal
mardbrdafullmat <- Anova(dbrdafullmat, type = 3)




#some script to print the names
namemar <- as.factor(row.names(subset(mardbrdafullmat, mardbrdafullmat[,3] <= 0.10)))

print(namemar)



VPAMF <- varpart(ddat, plantStructure, plantdiv, plantProd, soils)
VPAMF

plot(VPAMF)
windows()
bac1varpart <- 
  plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'plant productivity', 'Soils'), bg = c('navy', 'tomato','pink', 'goldenrod'))




#reduced marginal model based on sig value

marredmodmat<-dbrda(ddat ~ plantStructure + plantdiv  +     plantProd     + soils     + finalGEO  )

maraovredmodmat <- Anova(marredmodmat,type = 3)

# #models for individual variables

dbrdafullvar <-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)



marfullvar <- Anova(dbrdafullvar, type = 3)


#names for marginal model

namemar <- as.factor(row.names(subset(marfullvar, marfullvar[,3] <= 0.10)))

print(namemar)


#reduced marginal model based on sig value

marredmodvar<-dbrda(ddat ~ PlantStructure$vert.sd    +PlantStructure$VAI     +   PlantDiv$tree_richness  +  PlantDiv$AMdominance   +   PlantProd$stand_age_years + Soils$cn + Soils$pH     +     Soils$Fe_percent + finalGEO   )
maraovredmodvar <- Anova(marredmodvar,type = 3)


#names for sequential model
# 
# nameter <- as.factor(row.names(subset(terfullvar, terfullvar[,4] <= 0.050)))
# 
# print(nameter)
# 
# #reduced sequential model based on individual variables
# 
# terredmodvar<-dbrda(ddat ~ PlantStructure$vert.sd + PlantStructure$VCI   +  PlantDiv$AMdominance+   Soils$cn       +        Soils$pH         +      Soils$Fe_percent    + Condition(finalGEO))
# aovredmodter <- Anova(redmodter,by="terms", permu=999)
# 
maraovredmodvar$row.names <- row.names(maraovredmodvar)

maraovredmodmat$row.names <- row.names(maraovredmodmat)



#dataframe with all models
EMF1 <- data.frame(maraovredmodmat)

EMF1[0:7,5:9] <-   maraovredmodvar

#csv for easy copypasta
write.csv(EMF1, "EMF_0_5_chi.csv")







#-------------------- EMF 2----------------

setwd("/Users/joe/Desktop/Lang_etal_data_archive/data/betadiversity/EMF")
myfileOTU <- read.csv("IndianaCFIITSProp_seqtabEcto_5_10_2022_02_24.csv", head =TRUE, row.names=1)
myfileOTU<-t(myfileOTU)

PlantDiv <- read.csv("PlantTreeDiv_5_10_02_2022.csv", head =TRUE, row.names=1)
PlantStructure <- read.csv("PlantStructural_5_10_02_2022.csv", head =TRUE, row.names=1)
PlantProd <- read.csv("PlantProd_5_10_02_2022.csv", head =TRUE, row.names=1)
Soils <- read.csv("fe_pH_5_10_2022_02_15.csv", head =TRUE, row.names=1)
Geo <- read.csv("Geo_5_10_02_2022.csv", head =TRUE, row.names=1)
#myfileENV<-read.csv(file.choose(), head =TRUE, row.names=1)


ddat<-vegdist(as.matrix(myfileOTU), "jaccard")

detrend<-dbrda(ddat~., data=Geo)
anova(detrend,by="margin",permutations=how(nperm=9999))
resids<-resid(detrend)
mems.all<-dbmem(Geo, silent=FALSE)
memmod<-dbrda(resids~.,data=mems.all)
ordistep(memmod)
# this output will show you which MEMs to keep
#keeperMEMs<-c(mems.all[,c(4)])
keeperMEMs<-c(mems.all)

finalGEO<-as.matrix(cbind(Geo, keeperMEMs))
soils<-as.matrix(Soils)
plantdiv<-as.matrix(PlantDiv)
plantStructure<-as.matrix(PlantStructure)
plantProd<-as.matrix(PlantProd)
# 



dbrdafullmat <- dbrda(ddat ~ plantStructure + plantdiv + plantProd + soils + finalGEO)


#anova models

#marginal
mardbrdafullmat <- Anova(dbrdafullmat, type = 3)



#some script to print the names
namemar <- as.factor(row.names(subset(mardbrdafullmat, mardbrdafullmat[,3] <= 0.10)))

print(namemar)



VPAMF <- varpart(ddat, plantStructure, plantdiv, plantProd, soils)
VPAMF

plot(VPAMF)
windows()
bac1varpart <- 
  plot (VPAMF, digits = 3, Xnames = c('Plant Structure', 'Plant Diversity', 'plant productivity', 'Soils'), bg = c('navy', 'tomato','pink', 'goldenrod'))


#reduced marginal model based on sig value

marredmodmat<-dbrda(ddat ~ plantStructure + plantdiv      + plantProd +     soils  + finalGEO   )

maraovredmodmat <- Anova(marredmodmat,type = 3)
# 
# #models for individual variables

dbrdafullvar <-dbrda(ddat ~  PlantStructure$vert.sd + PlantStructure$VAI + PlantStructure$VCI + PlantDiv$tree_richness + PlantDiv$AMdominance + PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH + Soils$Fe_percent + finalGEO)


marfullvar <- Anova(dbrdafullvar, type = 3)
#names for marginal model

namemar <- as.factor(row.names(subset(marfullvar, marfullvar[,3] <= 0.10)))

print(namemar)


#reduced marginal model based on sig value

marredmodvar<-dbrda(ddat ~ PlantStructure$vert.sd  + PlantStructure$VAI     +   PlantDiv$tree_richness  +   +   PlantProd$stand_age_years + PlantProd$BAI  + Soils$cn + Soils$pH     +     Soils$Fe_percent + finalGEO   )
maraovredmodvar <- Anova(marredmodvar,type = 3)


maraovredmodvar$row.names <- row.names(maraovredmodvar)

maraovredmodmat$row.names <- row.names(maraovredmodmat)

#dataframe with all models
EMF2 <- data.frame(maraovredmodmat,  maraovredmodvar)


#csv for easy copypasta
write.csv(EMF2, "EMF_5_10_chi.csv")




