
#'
#' @title Get_Habitat_Data()
#' 
#' @description Get habitat data from Abell et al.'s (2008) freshwater ecoregion map
#' 
#' @details This function takes a data.frame with latitude and longitude data and gets
#' data on biogeographic realm, major habitat type and ecoregion based on Abell et al.'s (2008)
#' freshwater ecoregion map.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param data - data.frame with a column containing latitude and longitude in decimal degrees
#' @param latitude_dd - character string with the column name containing the latitude in decimal degrees
#' @param longitude_dd - character string with the column name containing the longitude in decimal degrees
#' 
#' @return data.frame with habitat data from Abell et al.s (2008) ecoregion map for a given set of coordinates
#' 

# function to get habitat data
Get_Habitat_Data <- function(data, latitude_dd, longitude_dd) {
  
  # load the freshwater habitat map
  if (!exists("fw_map")) {
    fw_map <- readRDS(file = here::here("database/freshwater_ecoregion_map.rds"))
  }
  
  # load the freshwater habitat map
  if (!exists("fw_meta")) {
    fw_meta <- readRDS(file = here::here("database/freshwater_ecoregion_metadata.rds"))
  }
  
  # make sure the correct packages are installed
  test_1 <- function(x) {
    
    all( c("dplyr", "sp", "raster") %in% installed.packages()[,1])
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(c("dplyr", "sp", "raster"), " must be installed for this function to work")
    
  }
  
  # test if the data input is a data.frame or a tibble
  test_2 <- function(x) {
    
    ( is.data.frame(x) | dplyr::is.tbl(x) )
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is not a data.frame or tibble object")
    
  }
  
  assertthat::assert_that(test_2(data))
  
  # test if the latitude and longitude columns are present in the data object
  test_3 <- function(x, y, z) {
    
    (assertthat::is.string(y) & assertthat::is.string(z)) & all(c(y, z) %in% names(x))
    
  }
  
  assertthat::on_failure(test_3) <- function(call, env){
    
    paste0(deparse(call$y), " or ", deparse(call$x),  " are not strings and/or are not present in the data object")
    
  }
  
  assertthat::assert_that(test_3(x = data, y = latitude_dd, z = longitude_dd))
  
  # test if the data object has more than 0 rows
  test_4 <- function(x) {
    
    (nrow(x) > 0)
    
  }
  
  assertthat::on_failure(test_4) <- function(call, env){
    
    paste0(deparse(call$x),  " object has zero rows and therefore no data")
    
  }
  
  assertthat::assert_that(test_4(x = data))
  
  # test if the data object has more than 0 rows
  test_5 <- function(x, y) {
    
    ( is.numeric(x) & is.numeric(y) ) | (all(is.na(x)) & all(is.na(x)))
    
  }
  
  assertthat::on_failure(test_5) <- function(call, env){
    
    paste0(deparse(call$x), " or ",deparse(call$y), " are not numeric variables")
    
  }
  
  assertthat::assert_that(test_5(x = data[[latitude_dd]], data[[longitude_dd]]))
  
  # make sure the decimal degree variables are within the ranges of the variables
  test_6 <- function(x, y) {
    
    t1 <- all( x[!is.na(x)] <= 90 )  
    t2 <- all( x[!is.na(x)] >= -90 )
    
    t3 <- all( y[!is.na(y)] <= 180 )
    t4 <- all( y[!is.na(y)] >= -180 )
    
    all( c(t1, t2, t3, t4) )
    
  }
  
  assertthat::on_failure(test_6) <- function(call, env){
    
    paste0(deparse(call$x), " or ",deparse(call$x), " are too big/small to be valid decimal degrees")
    
  }
  
  assertthat::assert_that(test_6(x = data[[latitude_dd]], data[[longitude_dd]]))
  
  # make sure the latitude-longitude data are numeric variables
  data[[latitude_dd]] <- as.numeric(data[[latitude_dd]])
  data[[longitude_dd]] <- as.numeric(data[[longitude_dd]])
  
  # add a row id column
  data$row_id <- 1:nrow(data)
  
  # remove the NA values
  lat.lon <- data[!(is.na(data[[longitude_dd]] ) | is.na(data[[latitude_dd]] )), ]
  
  # if there are no non-NA values then we don't get habitat data
  if( nrow(lat.lon) > 0 ) {
    
    # convert this to a spatial points object
    sp.pts <- sp::SpatialPoints(lat.lon[, c(longitude_dd, latitude_dd)], proj4string = raster::crs(fw_map) )
    
    # check where pts overlap with freshwater ecoregion
    hab.dat <- sp::over(sp.pts, fw_map)
    
    # add row id to these data
    hab.dat$row_id <- lat.lon$row_id
    names(hab.dat) <- c("habitat_id", "area_km2", "row_id")
    
    # join the habitat data to the original data
    data <- dplyr::full_join(data, hab.dat, by = "row_id")
    
  } else {
    
    data[["habitat_id"]] <- NA
    data[["area_km2"]] <- NA
    data[["row_id"]] <- NA
    
  }
  
  # arrange by row_id
  data <- dplyr::arrange(data, row_id)
  
  # remove the row_id column 
  data <- dplyr::select(data, -row_id, -area_km2)
  
  # join these data to the metadata
  data <- dplyr::left_join(data, fw_meta, by = "habitat_id")
  
  # convert to a tibble
  data <- dplyr::as_tibble(data)
  
  return(data)
  
}


