# for loading and wrangling rasters
library(raster)
# for transformations of shapefiles to raster crs
library(rgdal)
# for levelplot() function
library(rasterVis)
# for overlaying masked raster onto interactive basemap
library(leaflet)
# for saving leaflet maps
library(htmlwidgets)



# load map data as rasters
inf22 <- raster("data/Infra_Raster_2022-05-07.tiff", band = 1)
red22 <- raster("data/Infra_Raster_2022-05-07.tiff", band = 2)
inf23 <- raster("data/Infra_Raster_2023-05-12.tiff", band = 1)
red23 <- raster("data/Infra_Raster_2023-05-12.tiff", band = 2)
# the rasters 4 bands in total, though the last one is essentially empty
# only need first 2: infrared and red

# check crs
inf22@crs
inf23@crs

# creating raster from calculated NDVI, meaning: NDVI = (NIR - red) / (NIR + red)
NDVI_func <- function(inf_raster,red_raster) {
  NDVI <- (inf_raster - red_raster) / (inf_raster + red_raster)
  return(NDVI)
}

ndvi22 <- overlay(inf22, red22, fun = function(inf22, red22) NDVI_func(inf22, red22))
plot(ndvi22)
# oddly high values appear in the sea near harbor - possibly reflections or the like

# remove bodies of water by defining areas absorbing >90% of infrared as NA
ndvi22[inf22 < 20] <- NA

# same for 2023 .tiff
ndvi23 <- overlay(inf23, red23, fun = function(inf23, red23) NDVI_func(inf23, red23))
plot(ndvi23)

# less noise, but let's remove the water to be consistent
ndvi23[inf23 < 20] <- NA

# keeping only area covered by both maps
crop22 <- crop(ndvi22, extent(ndvi23))
crop23 <- crop(ndvi23, extent(ndvi22))

# set resolution to the smaller of the two files
crop22 <- resample(crop22, crop23, method = "bilinear")
crop23 <- resample(crop23, crop22, method = "bilinear")

# calculate differences between rasters from 2022 to 2023
rast_diff <- (crop23 - crop22) / 2
# we're only interested in NDVI changes in green areas
# so we filter out areas of NDVI less than 0 in both years
rast_diff[crop22 < 0 & crop23 < 0] <- NA



# creating a transparent color for mapping
transparent <- rgb(255, 255, 255, max = 255, alpha = 0)
# creating a color scale with transparency in the middle, making weaker effects fade from the map
col <- c("brown4","sandybrown", transparent, "green3", "green4")

levelplot(rast_diff, col.regions = colorRampPalette(col),
          at = seq(-1,1,length.out = 11),
          scales = list(draw = FALSE), margin = FALSE)
# the color balance is skewed by a handful of extreme changes

# these extremes can be capped to a smaller but still meaningful change
small_diff <- rast_diff
small_diff[small_diff > 0.2] <- 0.2
small_diff[small_diff < -0.2] <- -0.2

levelplot(small_diff, col.regions = colorRampPalette(col),
          at = seq(-0.2,0.2,length.out = 11),
          scales = list(draw = FALSE), margin = FALSE)
# the green shows where vegetation density grew going into 2023
# the brown shows where it shrank



# do the same with 2021, as that year was slightly wetter than 2023
# the satellite image is dated three weeks later than the 2022 one, however
inf21 <- raster("data/Infra_Raster_2021-05-29.tiff", band = 1)
red21 <- raster("data/Infra_Raster_2021-05-29.tiff", band = 2)

ndvi21 <- overlay(inf21, red21, fun = function(inf21, red21) NDVI_func(inf21, red21))
plot(ndvi21)
ndvi21[inf21 < 20] <- NA
  
crop22_2 <- crop(ndvi22, extent(ndvi21))
crop21 <- crop(ndvi21, extent(ndvi22))

crop22_2 <- resample(crop22_2, crop21, method = "bilinear")
crop21 <- resample(crop21, crop22_2, method = "bilinear")

levelplot(crop21, col.regions = colorRampPalette(col),
          at = seq(-1,1,length.out = 11),
          scales = list(draw = FALSE), margin = FALSE)
levelplot(crop22_2, col.regions = colorRampPalette(col),
          at = seq(-1,1,length.out = 11),
          scales = list(draw = FALSE), margin = FALSE)

# though the chronology is the same as above, this will map changes from '21->'22
rast_diff_2 <- (crop22_2 - crop21) / 2
rast_diff_2[crop21 < 0 & crop22_2 < 0] <- NA

levelplot(rast_diff_2, col.regions = colorRampPalette(col),
          at = seq(-1,1,length.out = 11),
          scales = list(draw = FALSE), margin = FALSE)

# making an interactive map out of it
m1 <- leaflet() %>%
  addTiles() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addRasterImage(rast_diff_2,
                 colors = col, 
                 opacity = 1)
m1
# here, the suburban areas again drown out the changes to the inner city

small_diff_2 <- rast_diff_2
small_diff_2[small_diff_2 > 0.2] <- 0.2
small_diff_2[small_diff_2 < -0.2] <- -0.2

levelplot(small_diff_2, col.regions = colorRampPalette(col),
          at = seq(-0.2,0.2,length.out = 11),
          scales = list(draw = FALSE), margin = FALSE)
# now, this shows how the NDVI changed from 2021 to 2022
# at a glance, while the drier year hurt much public land, gardens and the like seemed to thrive
# to examine this, we can focus on public green spaces: woods, parks, etc.

# first, making the raster less grainy before mapping
rast_dis_2 <- disaggregate(rast_diff_2, 2)
# smoothing out the local noise
rast_foc_2 <- focal(rast_dis_2, w=matrix(1/9, nrow=3, ncol=3))

# loading in geometry objects for public green spaces
parks <- readOGR("data/aarhus_parks/dynlayer5060129131777828764.shp")
forests <- readOGR("data/aarhus_forests/dynlayer3076536501278714013.shp")
yards <- readOGR("data/aarhus_churchyards/dynlayer10798847952077866102.shp")

# gathering the polygons in a single object for simpler masking
green_areas <- raster::union(parks, forests)
green_areas <- raster::union(green_areas, yards)

# reprojecting to raster CRS
# (the raster CRS is unfitting, but the raster itself becomes stretched when projected - leaflet will alleviate much of this)
green_proj <- spTransform(green_areas, crs(rast_dis_2))
# cropping the polygons to keep only those inside raster boundaries
green_crop <- crop(green_proj, rast_dis_2)

# creating mask for raster to show only forests, parks, and graveyards
rast_dry <- mask(rast_foc_2, green_crop)
# calibrating the maximum value intensities so individual values representing extreme highs/lows don't drown out the majority
rast_dry[rast_dry > 0.2] <- 0.2
rast_dry[rast_dry < -0.2] <- -0.2

levelplot(rast_dry, col.regions = colorRampPalette(col),
          at = seq(-0.2,0.2,length.out = 11),
          scales = list(draw = FALSE), margin = FALSE)

# creating interact-able map with leaflet
m2 <- leaflet() %>%
  addTiles() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addRasterImage(rast_dry,
                 colors = col, 
                 opacity = 1)
m2
# green-colored areas held more plant life in the dry 2022 than 2021 
# meanwhile, brown-colored areas held less plant life in the dry year 
# as such, the brown-colored areas are those that will suffer most in drought


# creating another map with inverted values to show where wet years create the most greenery
rast_diff_3 <- (crop21 - crop22_2) / 2
rast_dis_3 <- disaggregate(rast_diff_3, 2)
rast_foc_3 <- focal(rast_dis_3, w=matrix(1/9, nrow=3, ncol=3))
rast_wet <- mask(rast_foc_3, green_crop)

rast_wet[rast_wet > 0.2] <- 0.2
rast_wet[rast_wet < -0.2] <- -0.2


m3 <- leaflet() %>%
  addTiles() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addRasterImage(rast_wet,
                 colors = col, 
                 opacity = 1)
m3
# the green areas here are the ones most likely to be under threat by drought

# to check whether this effect is due to the later date at which the 2021 image was captured, wathe NDVI diff is mapped for 22-23
rast_dis <- disaggregate(small_diff, 2)
rast_foc <- focal(rast_dis, w=matrix(1/9, nrow=3, ncol=3))

green_proj_2 <- spTransform(green_areas, crs(rast_dis))
green_crop_2 <- crop(green_proj_2, rast_dis)
rast_23 <- mask(rast_foc, green_crop)

m4 <- leaflet() %>%
  addTiles() %>%
  addProviderTiles(providers$CartoDB.Positron) %>% 
  addRasterImage(rast_23,
                 colors = col, 
                 opacity = 1)
m4

saveWidget(m2, file="out/map_2021_to_22.html")
saveWidget(m3, file="out/map_2022_to_21.html")
saveWidget(m4, file="out/map_2022_to_23.html")
