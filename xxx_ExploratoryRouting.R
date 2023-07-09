

library(tidyverse)
library(sf)
library(sfnetworks)




library(marmap)
library(raster)
library(sf)

#Make three ship positions near PNG
pts <- data.frame(nr=1:3, x=c(145,145,150), y=c(-10,-1,-7)) %>% 
  sf::st_as_sf(coords = c("x","y")) %>% 
  sf::st_set_crs(4326)

#Get bathymetric data from the area
papoue <- getNOAA.bathy(lon1 = 140, lon2 = 155,
                        lat1 = -13, lat2 = 0, resolution = 4)

#Turn into raster and keep only depth lower than 0 meter (water)
r <- marmap::as.raster(papoue)
plot(r)

r@data@values[r@data@values>=0] <- NA
r@data@values[r@data@values<0] <- 1

plot(r)
plot(pts, add=T)