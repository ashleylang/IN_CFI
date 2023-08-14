#IN CFI stand-level attributes
#Jan. 2021

#data
wd <- "/data/"
tree <- read.csv(paste(wd, "/CFI_inventory.csv", sep=""))

#Calculate stand attributes from CFI tree level data that match plot names in final dataset
#Plot coords and tree level data using plotID from 2015 and 2020 inventories must be requested from Indiana DNR

#Trees > 5 inch or 12.7 cm DBH - which is the minimum size that is sampled across the whole plot
#calculate diameter in cm and basal area of every tree
#m2 = (dbh in cm ^2) * 0.00007854 
tree$DIA_CM2020 <- tree$diameter2020_in * 2.54
tree$DIA_CM2015 <- tree$diameterprevious_in * 2.54
tree$basal_area2020 <- (tree$DIA_CM2020^2) * 0.00007854
tree$basal_area2015 <- (tree$DIA_CM2015^2) * 0.00007854

#############################################################################################################
#Calculate stand attributes with a loop to zip all the data together quickly (only ~ 40 plots)
#############################################################################################################
#old fashioned loop to calculate stand attributes for each plot
plots <- unique(tree$county_plotID)

OUT <- NULL
for(n in 1:length(plots)) {
  plot.name.i <- plots[n]
  tree.i <- tree[tree$county_plotID == plot.name.i, ]
  #2020 data (must remove dead trees b/c 2015 data is nested within dataset)
  tree.i.2020 <- tree.i[tree.i$status2020 == "live", ]
  
  #tree species richness in 2020
  tree_richness <- length(unique(tree.i.2020$SPCD))
  
  #total basal area
  BAtotal_2020 <- sum(tree.i.2020$basal_area2020, na.rm = TRUE)
  
  #basal area increments
  #2015 data (must only have live 2015 trees)
  tree.i.2015 <- tree.i[tree.i$statusprevious == "live", ]
  BAtotal_2015 <- sum(tree.i.2015$basal_area2015, na.rm = TRUE)
  BAI <- (BAtotal_2020 - BAtotal_2015) / 5
  
  #EM and AM tree species richness and basal area
  EM <- tree.i.2020[tree.i.2020$EM == 1, ]
  EMtree_richness <- length(unique(EM$SPCD))
  EMBA <- sum(EM$basal_area2020, na.rm = TRUE)
  
  AM <- tree.i.2020[tree.i.2020$AM == 1, ]
  AMtree_richness <- length(unique(AM$SPCD))
  AMBA <- sum(AM$basal_area2020, na.rm = TRUE)
  
  AMdominance <- AMBA / (AMBA + EMBA)
  
  out.i <- cbind.data.frame(plot.name.i, tree.i[1,2:7], tree_richness,
                            BAtotal_2020, BAI, EMtree_richness, EMBA, AMtree_richness, 
                            AMBA, AMdominance)
  OUT <- rbind(OUT, out.i)
  
}

write.csv(OUT, paste(wd, "./CFI_soil_sampling_stand_attributes_01132021.csv", sep=""))