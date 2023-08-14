################################################################################
#Obtain USGS 3DEP lidar, clean it, and calculated structural diveristy metrics
#Feb. 2021
#this works with the Version 3.0.3 of lidR (8/3/2020)
################################################################################
setwd("./data/")

######################## Packages ##############################################
library('lidR') 
#####################################################################################
#can be freely obtained from USGS 3DEP Program around study coordinates provided
ctgE <- readLAScatalog(folder = "./variables/lidar/INCFIPlots_1500/eastCFI/")
ctgW <- readLAScatalog(folder = "./variables/lidar/INCFIPlots_1500/westCFI/")
#####################################################################################
# Part of this must be done in ArcGIS beforehand. Indiana uses a unique coordinate system ##
## that divides the state into east/west halves that don't overlap. You cannot change the CRS for LiDAR ##
## files easily, so you must instead change the shapefile containing plot coordinates to the same CRS as the ##
## LiDAR data. ##

## Create an excel file with x and y coordinates for the plot centers in ArcMap
## one from the Indiana state plane east projection and one from west. ##
shp.east <- shapefile("./variables/lidar/GIS/myCFI_INEast/myCFI_INEast.shp")
shp.west <- shapefile("./variables/lidar/GIS/myCFI_InWest/myCFI_InWest.shp") 

##The initial clip from the USGS database created a 1500 ft buffer around each plot. 
#Note lidar spatial units in ft
#here we clip down to 100 ft to make processing more manageable
opt_output_files(ctgE) <- "./variables/lidar/Processed/eastCFI/eastCFI_buffer/{county_plo}"
east <- lasclip(ctgE, shp.east, radius = 100) 

opt_output_files(ctgW) <- "./variables/lidar/Processed/westCFI/westCFI_buffer/{county_plo}"
west <- lasclip(ctgW, shp.west, radius = 100)

#####################################################################################
#####################################################################################
#Remove outliers and correct point height to account for ground elevation ##
#####################################################################################
#function calculate to and remove unwanted outliers in readLAS()
outliers <- function(las_file) {
  #remove outlier points using mean and sd statistics
  Z_sd=sd(las_file@data$Z)
  Z_mean=mean(las_file@data$Z)
  #using 6 sd as a coarse filter due to topopographic variation
  #will fine filter later after corrected for elevation if visually outliers are found
  f= paste("-drop_z_below",(Z_mean-6*Z_sd),"-drop_z_above",(Z_mean+6*Z_sd))
  return(f)
}
# --------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#indiana lidar uses Indiana East and West state plane so plots split 
#indicate which side processing and then run loop for each side
side <- "westCFI"
side <- "eastCFI"

dir.create(paste("./variables/lidar/Processed/", side, "/rm_noise/", sep=""))
buffer.names <- list.files(paste("./variables/lidar/Processed/", side, "/", side, "_buffer/", sep="")) 

for(i in 1:length(buffer.names)){
  #reading in individual las plots not working with a catalog here
  las_file<- readLAS(file.path(paste("./variables/lidar/Processed/", side, "/", side, "_buffer/", buffer.names[i], sep="")))
  print(buffer.names[i])
  print(range(las_file@data$Z))
  
  #drop out outliers and rewrite a new file
  f <- outliers(las_file)
  las_file <- readLAS(file.path(paste("./variables/lidar/Processed/", side, "/", side, "_buffer/", buffer.names[i], sep="")), filter = f)
  print(range(las_file@data$Z))
  #correct for ground height
  # note function removes deprecated points which is indicated in that this happens in error message
  las_file <- normalize_height(las_file, tin())
  
  writeLAS(las_file, paste("./variables/lidar/Processed/", side, "/rm_noise/", buffer.names[i], sep=""))
}

# --------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
rm_noiseE <- readLAScatalog(folder = "./variables/lidar/Processed/eastCFI/rm_noise/")
las_check(rm_noiseE)
plot(rm_noiseE)
                                                  
rm_noiseW <- readLAScatalog(folder = "./variables/lidar/Processed/westCFI/rm_noise/")
las_check(rm_noiseW)
plot(rm_noiseW)
##################################################################################### 
#####################################################################################
#Clip final 24 ft radius plot extent to match CFI tree plot radius
#####################################################################################
#####################################################################################
#This function clips down to final 7.3 m  or 24 ft radius for structural metrics. 
#Make sure to use same shapefile as the initial clip above

opt_output_files(rm_noiseE) <- "./variables/lidar/processed/eastCFI/east24ftradius/{county_plo}"
east <- lasclip(rm_noiseE, shp.east, radius = 24) 

opt_output_files(rm_noiseW) <- "./variables/lidar/processed/westCFI/west24ftradius/{county_plo}"
west <- lasclip(rm_noiseW, shp.west, radius = 24)
#manually checked each file for outliers 

#####################################################################################
#####################################################################################
#Forest structural diversity function
#converted from ft to m in the function
############################################################

structural_diversity_metrics <- function(data.40m) {
  
  #metrics from cloud_metrics (convert ft to m post calc)
  vert.sd <- cloud_metrics(data.40m, sd(Z, na.rm = TRUE)) 
  vert.sd <- vert.sd * .3048
  
  #metrics from a Z vector
  Zs <- data.40m@data$Z
  Zs <- Zs[!is.na(Zs)]
  Zs <- Zs * 0.3048 #just convert Z vector from ft to m
  LADen<-LAD(Zs, dz = 3.28, k=.5, z0=9.84) 
  VAI <- sum(LADen$lad, na.rm=TRUE) 
  VCI <- VCI(Zs, by = 3.28, zmax=200) 
  out.plot <- data.frame(matrix(c(vert.sd, VAI,
                                  VCI), ncol = 3)) 
  colnames(out.plot) <- c("vert.sd",
                          "VAI", "VCI")
  return(out.plot)
}


#####################################################################################
#loop to generate metric table 
#################################################################################
OUT <- NULL #only once and then combine e and w into one table

#select the side of the state and run through each once
las.names <- list.files("./variables/lidar/Processed/eastCFI/east24ftradius/")
las.folder <- "./variables/lidar/Processed/eastCFI/east24ftradius/"
las.names <- list.files("./variables/lidar/Processed/westCFI/west24ftradius/") 
las.folder <- "./variables/lidar/Processed/westCFI/west24ftradius/"

#i=1
for(i in 1:length(las.names)){
  #-----------------------------------------------------------------
  #before filtering ground points for quality check ----------------
  #reading in individual las plots without filtering ground points
  data.40m <- readLAS(file.path(paste(las.folder, las.names[i], sep="")))
  print(las.names[i])
  
  #get area so can filter out cropped plots #convert ft2 to m2
  plot_area_ground <- area(data.40m) * 0.092903
  
  #total number of points divided by plot area
  den <- length(data.40m@data$Z)/plot_area_ground
  
  #------------------------------------------------------------------
  #reading in individual las plots again with filtering outliers and points < 2 m
  data.40m <- readLAS(file.path(paste(las.folder, las.names[i], sep="")), filter = "-drop_z_below 6.561")
  
  #remove duplicated ("deprecated points") for heterogeneity metrics
  data.40m <- filter_duplicates(data.40m)
  
  #get area of veg points #convert ft2 to m2
  veg_plot_area <- area(data.40m) * 0.092903
  
  #if there are no points because all veg is < 2 m keep loop running with if else options
  if(veg_plot_area > 0) { #if veg > 2 m 
    #run the FSD metric function
    FSD.i <- structural_diversity_metrics(data.40m)
    
  } else { #if veg < 2 m
    FSD.i <- data.frame(matrix(rep(0, 3), ncol = 3))
    colnames(FSD.i) <- c("vert.sd",
                         "VAI", "VCI")
  }
  
  #-------------------------------------------------------------------
  out.i <- cbind.data.frame(las.names[i], plot_area_ground, 
                            den, FSD.i)
  OUT <- rbind(OUT, out.i)
}

colnames(OUT)[1] <- "plotID"

write.csv(OUT, file = "./variables/CFI_myc_lidarFSD_05042021_24ftradius.csv")

#####################################################################################

#End Workflow

#####################################################################################

