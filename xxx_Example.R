



# User input --------------------------------------------------------------


precip_threshold <- 1.5 #mm hr-1
storm_radius <- 250 # threshold in km

# End user input ----------------------------------------------------------

library(sf)
library(terra)
library(tidyverse)

source('read_track.r')
source('associate_storm.r')


# Load sample data --------------------------------------------------------

storms <- st_read("data/EXAMPLE_track_timestep.geojson")
interest <- st_read("data/england_wales_box.geojson")
precip <- terra::rast("data/EXAMPLE_era5_precip.tif")


# Run association function ------------------------------------------------

assoc_results <- associate_storm(storms, precip, interest, precip_threshold, storm_radius)




