#Create a site map
#March 2022

#############################################################################
#the data ---------------------------------------------------------------------
setwd("./data/")
dat <- read.csv("./Geo_0_5_02_2022.csv")
dat <- dat[, c("GPS_LON", "GPS_LAT")] #keep latitude and longitude columns
dat <- dat[!duplicated(dat),] #remove any duplicated points
##############################################################################
library(tmap) 
library(sf)
library(USAboundaries)
library(ggplot2)
library(cowplot)
#############################################################################
#Figure 2 ------------------------------------------------------
#Indiana boundary
ind <- us_states() #call the us state borders from USAboundaries package
remove <- c("AK", "HI", "PR")
ind <-ind[!(ind$stusps %in% remove),]

ind2 <- ind[ind$stusps == "IN", ]
st_crs(ind)

#create sf object
site_sf <- st_as_sf(dat, coords = c("GPS_LON", "GPS_LAT"), crs = st_crs(ind)) 

#transform to a nicer looking projection
ind <- st_transform(ind, crs = st_crs(5070)) 
ind2 <- st_transform(ind2, crs = st_crs(32616))
site_sf <- st_transform(site_sf, crs = st_crs(32616)) 

#two panel map
tm1 <- tm_shape(ind) + tm_borders("black", lwd = .5) + tm_shape(site_sf) + 
  tm_dots(size = .1, shape = 3, col = "black") + tm_scale_bar()
tm2 <- tm_shape(ind2) + tm_borders("black", lwd = .5) + tm_shape(site_sf) + 
  tm_dots(size = .1, shape = 4, col = "black") + tm_layout(
    inner.margins = c(.3, .3, 0, 0), frame = FALSE)
tmap_arrange(tm1, tm2, ncol = 2, nrow = 1)


#inset map
main_tmap <- tmap_grob(tm2)
sub_tmap <- tm_shape(ind) + tm_polygons() + tm_layout(frame = FALSE)
sub_tmap <- tmap_grob(sub_tmap)

ggdraw() +
  draw_plot(main_tmap) +
  draw_plot(sub_tmap,
            height = 0.3,
            x = 0.05)
