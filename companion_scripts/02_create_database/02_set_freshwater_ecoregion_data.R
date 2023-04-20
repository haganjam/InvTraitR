#' @title freshwater ecoregion database
#' @description create a data.frame of realms, habitat types and ecoregions
#' @details This script uses the freshwater ecoregion map (Abell et al. 2008)
#'  to get realm, major habitat type and eco-region data for each equation,
#'  length
#' etc. in the database.
#' @author James G. Hagan (james_hagan(at)outlook.com)

# load relevant spatial libraries
library(sp)
library(sf)
library(raster)
library(dplyr)

# load the set of metadata associated with each FEOW_ID
fw.md <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver4/freshwater_ecoregion_habitat_list.xlsx")
head(fw.md)

# remove the page column
fw.md <-
  fw.md %>%
  dplyr::select(-page)

# save the habitat map as .rds file as well
saveRDS(fw.md, file = paste("database", "/", "freshwater_ecoregion_metadata.rds", sep = ""))

# set-up the CRS
crdref <- CRS("+proj=longlat +datum=WGS84")

# load the freshwater map
fw <- st_read("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver4/feow_hydrosheds.shp")

# set the crs
st_crs(fw) <- crdref

# convert to spatial object
fw <- as(fw, "Spatial")

# load the habitat data for each equation, length etc.
hab <- readxl::read_xlsx(path = "C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_4_FreshInvTraitR/data/allometry_database_ver4/habitat_database.xlsx")
head(hab)

# create a spatial points object for each lat lon in the habitats database
df <- data.frame(
  longitude = as.numeric(hab$lon_dd),
  latitude = as.numeric(hab$lat_dd)
)
head(df)
str(df)
summary(df)

# add a row id column
df$row_id <- 1:nrow(df)

# remove the NA values
df2 <- df[!is.na(df$longitude), ]

# convert this to a spatial points object
pts <- SpatialPoints(df2[, c(1, 2)], proj4string = crs(fw))
plot(pts)

# check where pts overlap with freshwater ecoregion
my_over_output <- over(pts, fw)

# add row id to these data
my_over_output$row_id <- df2$row_id
names(my_over_output) <- c("habitat_id", "area_km2", "row_id")

# join to the original df to refill in the NAs
my_over_output <-
  full_join(my_over_output, df[3], by = "row_id") %>%
  arrange(row_id)

# remove the row_id column
my_over_output <-
  my_over_output %>%
  dplyr::select(-row_id)

# join these data to the metadata
my_over_output <- left_join(my_over_output, fw.md, by = "habitat_id")

# add the habitat data to the habitat database
hab$habitat_id <- my_over_output$habitat_id

# add realm data
hab$realm <- my_over_output$realm

# add major_habitat_type
hab$major_habitat_type <- my_over_output$major_habitat_type

# add ecoregion
hab$ecoregion <- my_over_output$ecoregion

# add area
hab$area_km2 <- my_over_output$area_km2

# reorganise the and arrange the columns
hab <-
  hab %>%
  dplyr::select(
    database, id, accuracy, lat_dd, lon_dd, habitat_id, area_km2,
    realm, major_habitat_type, ecoregion
  ) %>%
  arrange(database, id)
head(hab)

# write this into an rds file
saveRDS(hab, file = paste("database", "/", "freshwater_ecoregion_data.rds", sep = ""))

# save the habitat map as .rds file as well
saveRDS(fw, file = paste("database", "/", "freshwater_ecoregion_map.rds", sep = ""))

