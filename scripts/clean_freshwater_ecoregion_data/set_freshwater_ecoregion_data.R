
# incorporate the Freshwater ecoregions map

# load relevant spatial libraries
library(sp)
library(sf)
library(raster)

# set-up the CRS
crdref <- CRS('+proj=longlat +datum=WGS84')

# load the freshwater map
fw <- st_read("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_BEF_rockpools_Australia/data/trait_and_allometry_data/freshwater_ecoregion_map/GIS_hs_snapped/feow_hydrosheds.shp")

# set the crs
st_crs(fw) <- crdref

# convert to spatial object
fw = as(fw, 'Spatial')

# set-up a spatial points polygon

# define some latitude and longitudes
longitude <- c(-116.7, -120.4, -116.7, -113.5, -115.5,
               -50, 20, -100, -80, 10)
latitude <- c(45.3, 42.6, 38.9, 42.1, 35.7, 38.9,
              0, 60, 10, 50)

df <- data.frame(longitude, latitude)
head(df)

# convert this to a spatial points object
pts <- SpatialPoints(df, proj4string = crs(fw) )
pts
plot(pts)

# check where pts overlap with freshwater ecoregion
my_over_output = over(pts[1:5,], fw)
my_over_output

### END
