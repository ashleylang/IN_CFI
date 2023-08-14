############# extract MAT and MAP from CHELSA climate data for CFI plots
#Feb. 2021

#####################
library(raster)
library(sp)
#####################
setwd("./data/")
plots <- read.csv("./coordinates.csv")
options(stringsAsFactors = FALSE)
#Data downloaded from https://chelsa-climate.org - CHELSA Version 1-2
bio01 <- raster("./chelsa_02182021/CHELSA_bio10_01.tif")
bio12 <- raster("./chelsa_02182021/CHELSA_bio10_12.tif")


coords <- data.frame(x=plots$GPS_LON, y=plots$GPS_LAT)
rownames(coords) <- plots$county_plotID
points <- SpatialPoints(coords, proj4string = bio01@crs)


plot_bio01 <- extract(bio01, points)
plot_bio01 <- plot_bio01 / 10
plot_bio12 <- extract(bio12, points)

df <- cbind.data.frame(coordinates(points), plot_bio01, plot_bio12)
df$county_plotID <- rownames(df)

write.csv(df, "./climate_CHELSA_INCFI_plots_02192021.csv")
