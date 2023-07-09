

library(raster)
library(tidyverse)
library(sf)
library(mapiso)
library(gdistance)
library(geosphere)
library(ncdf4)

source('read_track.r')



# User input --------------------------------------------------------------


target_datetime <- "2014.02.04.18.00.00"
precip_threshold <- 1.5 #mm hr-1
track_threshold <- 250 # threshold in km
interest_area <- data.frame(xmin = -8, xmax = 3, ymin = 49, ymax = 55)
track_file <- "C:/Users/ruari/OneDrive/Documents/R_scripts/Data/Tracks/ERA-5/DJF/NH_FILT/ERA5_20132014_TRACKS_FILTERED_pos"
track_season_start <- "2013-12-01"
# End user input




# Load track file ---------------------------------------------------------

track_df <- track_file %>% 
  read_track() %>% 
  add_timestamp_to_tracks(start_date = "2013-12-01") %>% 
  convert_track_to_negative_longitude() %>% 
  crop_tracks_to_extent(extent = c(xmin = -70, xmax = 25, ymin = 20, ymax = 75))

track <- track_to_sf_points(track_df)
track_lines <- track_to_sf_lines(track_df)


# Find track timestep
track_timestep <- filter(track, time == "18:00", date == "2014-02-04")




# Interest area -----------------------------------------------------------

# Create polygon for interest_area
interest_area_df <- data.frame(
  x = with(interest_area, c(xmin, xmin, xmax, xmax, xmin)),
  y = with(interest_area, c(ymin, ymax, ymax, ymin, ymin))
)

interest_area_wkt <- sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))",
                             interest_area_df[1,1],
                             interest_area_df[1,2],
                             interest_area_df[2,1],
                             interest_area_df[2,2],
                             interest_area_df[3,1],
                             interest_area_df[3,2],
                             interest_area_df[4,1],
                             interest_area_df[4,2],
                             interest_area_df[5,1],
                             interest_area_df[5,2])
interest_area <- st_as_sf(data.frame(id = 1, geometry = interest_area_wkt), wkt = "geometry", crs = 4326)

interest_centroid <- st_centroid(interest_area)




# Precip data -------------------------------------------------------------

# Load precip data (if not yet loaded)
if(!exists("data_loaded")) data_loaded <- FALSE

target_datetime <- paste0("X", target_datetime)

# Load data only if required
if(!data_loaded){
  era5_data <- brick("C:/Users/ruari/OneDrive/Documents/R_scripts/Data/Reanalysis/ERA-5/era5-precip-2014-janfeb.nc")
  data_loaded <- TRUE
}

# Filter precip data to target date/time
era5_data_slice <- era5_data[[target_datetime]]


# Coerce raster brick to points
era5_pts <- as(era5_data_slice, "SpatialPixelsDataFrame")
era5_df <- as.data.frame(era5_pts) %>% 
  pivot_longer(-c("x","y"), names_to = "timestamp", values_to = "precip")


# units are kg m-2 s-1 (mm s-1)
# convert to mm/hour
era5_df$precip <- era5_df$precip * 3600

# Create 10mm/day contour
precip_contour <- mapiso(era5_df, var = "precip", coords = c("x","y"), breaks = precip_threshold, crs = "EPSG:4326")
# Exclude the < 10mm/day contour
precip_contour <- subset(precip_contour, isomin == precip_threshold)


precip_brick <- as(era5_data, "SpatRaster")

ext(precip_brick)



# Base map ----------------------------------------------------------------
# get world map data

library(tidyverse)
world_map <- map_data("world")

# Crop to extent
field_extent <- c(ext(precip_brick)[1], ext(precip_brick)[2],
                  ext(precip_brick)[3], ext(precip_brick)[4])
base_map <- filter(world_map, 
                long >= field_extent["xmin"], 
                long <= field_extent["xmax"], 
                lat >= field_extent["ymin"], 
                lat <= field_extent["ymax"])

# Remove polys with less than 3 entries
exclude_polys <- base_map %>% 
  group_by(group) %>% 
  tally() %>% 
  filter(n < 3) %>% 
  pull(group)

base_map <- filter(base_map, !(group %in% exclude_polys))



base_map_sf <- sfheaders::sf_polygon(
  obj = base_map,
  x = "long",
  y = "lat",
  polygon_id = "group"
)

# Remove polys with less than 3 entries
exclude_polys <- world_map %>% 
  group_by(group) %>% 
  tally() %>% 
  filter(n < 3) %>% 
  pull(group)

world_map <- filter(world_map, !(group %in% exclude_polys))

world_map_sf <- sfheaders::sf_polygon(
  obj = world_map,
  x = "long",
  y = "lat",
  polygon_id = "group"
)
base_map_2 <- st_crop(world_map_sf, ext(precip_brick))
st_crs(base_map_2) <- 4326

ggplot(base_map_2) + geom_sf()


ggplot() + geom_map(data = world_map_s, map = base_map,
                    aes(long, lat, map_id = region))




# Stormtrack data ---------------------------------------------------------


# Crop stormtracks to precip data region ----------------------------------

storm_tracks <- st_crop(track_timestep, extent(era5_data))

# Create 250km buffer around tracks
storm_tracks_lambert <- st_transform(storm_tracks, "EPSG:3035")
storm_tracks_buffer <- st_buffer(storm_tracks_lambert, track_threshold * 1000)

# Reproject buffer back to lat/lon for consistency
storm_tracks_buffer <- st_transform(storm_tracks_buffer, "EPSG:4326")







# Create a mask layer describing location of fronts, storm tracks, --------
# and target region

# Create mask layer to generate paths through fronts and track buffers
mask_base <- era5_data_slice
mask_base[,] <- 1

# Add precip contours to mask
mask_precip <- mask(mask_base, precip_contour)
mask_precip[is.na(mask_precip)] <- 0

# Add track buffers to mask
mask_track <- mask(mask_base, storm_tracks_buffer)
mask_track[is.na(mask_track)] <- 0

# Add interest area to mask
mask_interest <- mask(mask_base, interest_area)
mask_interest[is.na(mask_interest)] <- 0

# Combine mask layers
mask_full <- mask_precip + mask_track + mask_interest
mask_full[mask_full > 0] <- 1
mask_full[mask_full == 0] <- NA

# Find nearest connected storm centre -------------------------------------

# Create SpatialPoints object with point locations
# Need to automate this out of interest_centroid and storm_tracks


# Convert storm tracks to SPDF for compatibility with gDistance
track_pts <- as(storm_tracks, "Spatial")
interest_pt <- as(interest_centroid, "Spatial")


# Create a transition layer using gdistance package
tl <- transition(mask_full, transitionFunction = mean, directions = 4)

# Attempt to find the shortest valid path from the interest region centroid
# to each storm track point. If no valid path can be found, NA is returned
paths <- list()
for(i in 1:nrow(track_pts)){
  paths[[i]] <- tryCatch(
    shortestPath(tl, interest_pt, track_pts[i,], output = "SpatialLines"),
    error = function(e){return(NA)})
}

candidate_pts <- track_pts[which(!is.na(paths)),]
candidate_paths <- paths[!is.na(paths)]

# If multiple possible tracks are found, get their lengths and return the shortest
line_lengths <- list()
if(length(candidate_paths) > 1){
  for(i in 1:length(candidate_paths)){
    line_lengths[[i]] <- lengthLine(candidate_paths[[i]])
  } 
} else {
  line_lengths[[1]] <- lengthLine(candidate_paths[[1]])
}
line_lengths <- unlist(line_lengths)
shortest_line <- which.min(line_lengths)
matched_pt <- candidate_pts[shortest_line,]
matched_path <- st_as_sf(candidate_paths[[shortest_line]])



# Draw map -------------------------------------------------------------
era5_df %>% 
ggplot() +
  geom_tile(aes(x=x, y=y, fill =  precip)) +
  geom_map(
    data = world, map = world,
    aes(map_id = region),
    fill = NA, colour = "grey80"
  ) +
  geom_sf(data = precip_contour, colour = "red", fill = NA) + 
  geom_sf(data = storm_tracks, colour = "orange", size = 3) + 
  geom_sf(data = storm_tracks_buffer, colour = "orange", size = 2, fill = NA) + 
  geom_sf(data = interest_area, colour = "purple", size = 2, fill = NA) +
  geom_sf(data = matched_path) +
  scale_fill_viridis_c("Total precip (mm/hour)")#, limits = c(0, 1000))

