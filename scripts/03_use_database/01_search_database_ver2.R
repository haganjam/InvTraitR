
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
#' @return tibble of the input data with habitat data from Abell et al.s (2008) ecoregion map attached as additional columns
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

#'
#' @title Clean_Taxon_Names()
#' 
#' @description Get taxonomic distances of target names relative to the taxa databases
#' 
#' @details This function takes a data.frame with target_taxa and uses the bdc package
#' to clean and harmonise the names.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param data - data.frame with a column containing target taxon names and life-stages
#' @param target_taxon - character string with the column name containing the taxon names
#' @param life_stage - character string with the column name containing the life_stages
#' @param database - taxonomic database to use: gbif, itis, col
#' 
#' @return data.frame with target taxon names and cleaned names for the chosen taxonomic backbone
#' 

Clean_Taxon_Names <- function(data, target_taxon, life_stage, database = "gbif") {
  
  # set-up a vector of required packages
  pack_vec <- c("bdc", "dplyr", "taxadb", "curl")
  
  # make sure the correct packages are installed
  test_1 <- function(x) {
    
    all(pack_vec %in% installed.packages()[,1])
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(pack_vec, " must be installed for this function to work")
    
  }
  
  assertthat::assert_that(test_1())
  
  # test if the database input is a supported taxonomic backbone
  test_2 <- function(x) {
    
    assertthat::are_equal(x, "gbif") | 
      assertthat::are_equal(x, "itis") | 
      assertthat::are_equal(x, "col")
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is not a valid taxonomic backbone, pick: gbif, itis or col")
    
  }
  
  assertthat::assert_that(test_2(database))
  
  # test if the data input is a data.frame or a tibble
  test_3 <- function(x) {
    
    ( is.data.frame(x) | dplyr::is.tbl(x) )
    
  }
  
  assertthat::on_failure(test_3) <- function(call, env){
    
    paste0(deparse(call$x), " is not a data.frame or tibble object")
    
  }
  
  assertthat::assert_that(test_3(data))
  
  # test if the target_taxon column is in the data object
  test_4 <- function(x, y) {
    
    assertthat::is.string(y) & ( y %in% names(x) )
    
  }
  
  assertthat::on_failure(test_4) <- function(call, env){
    
    paste0(deparse(call$y), " is not a column in the supplied data object")
    
  }
  
  assertthat::assert_that(test_4(x = data, y = target_taxon))
  
  # test if the name column has a length of more than zero and that it is a character variable
  test_5 <- function(x) {
    
    is.character(x) & (length(x) > 0 )
    
  }
  
  assertthat::on_failure(test_5) <- function(call, env){
    
    paste0(deparse(call$x), " is not a character variable with length greater than zero")
    
  }
  
  assertthat::assert_that(test_5(data[[target_taxon]]))
  
  # test if the life-stage column is a character vector without NAs
  test_6 <- function(x) {
    
    (is.character(x) & all(x %in% c(NA, "none", "larva", "pupa", "nymph", "adult", "nauplius", "copepodite", "tadpole")))
    
  }
  
  assertthat::on_failure(test_6) <- function(call, env){
    
    paste0(deparse(call$x), " one or more entries do not have appropriate life-stage classes: see documentation")
    
  }
  
  assertthat::assert_that(test_6(data[[life_stage]]))
  
  # add a targ_no column to the data input
  data[["row_id"]] <- 1:nrow(data)
  
  # clean the names for typos etc. using the bdc_clean_names function
  clean.names <- bdc::bdc_clean_names(sci_names = data[[target_taxon]], save_outputs = FALSE)
  
  # write some code to remove the output file
  unlink("Output", recursive=TRUE)
  
  # add the clean names to the data.frame
  clean.col <- paste("clean_", target_taxon, sep = "")
  data[[clean.col]] <- clean.names$names_clean
  
  # check that all the non-missing names were cleaned
  test_7 <- function(x, y) {
    
    length(x) == length(y)
    
  }
  
  assertthat::on_failure(test_7) <- function(call, env){
    
    "Length of clean names do not match length of original names"
    
  }
  
  assertthat::assert_that(test_7(x = data[[target_taxon]], y = data[[clean.col]]))
  
  # update the database if there is a valid internet connection
  if (curl::has_internet()) {
    
    taxadb::td_create(provider = database,
                      overwrite = FALSE
    )
    
  }
  
  # add the database to the data
  data$db <- ifelse(is.na(data[[clean.col]]), NA, database )
  
  # subset out taxa with special names
  source(here::here("scripts/01_special_names_func.R"))
  spec.names <- special_taxon_names()
  data.spec <- dplyr::filter(data, (eval(parse(text = clean.col)) %in% spec.names) )
  
  # change the database column to special
  if (nrow(data.spec) > 0) {
    data.spec[["db"]] <- "special"
  }
  
  # remove the special names from the data
  data <- dplyr::filter(data, !(row_id %in% data.spec[["row_id"]]) )
  
  # if the are data points that are not special names, then we clean those names
  if (nrow(data) > 0) {
    
    # harmonise the names to the chosen data.base
    data.harm <- 
      bdc::bdc_query_names_taxadb(sci_name = data[[clean.col]],
                                  db = database,
                                  rank_name = "Animalia",
                                  rank = "kingdom"
      )
    
    # write some code to remove the output file
    unlink("Output", recursive=TRUE)
    
    # add a row_id to this harm.tax object
    data.harm$row_id <- data$row_id
    
    # higher taxon rank
    data.harm <- dplyr::mutate(data.harm,
                              db_taxon_higher_rank = ifelse(is.na(order) & is.na(family), NA, 
                                                            ifelse(is.na(order) & !is.na(family), "family", "order") ) )
    # higher taxon name
    data.harm <- dplyr::mutate(data.harm,
                              db_taxon_higher = ifelse(is.na(order) & is.na(family), NA, 
                                                       ifelse(is.na(order), family, order) ) )                          
    
    # select the relevant columns
    data.harm <- data.harm[,c("row_id", "scientificName", "acceptedNameUsageID", "db_taxon_higher_rank", "db_taxon_higher")]
    
    # remove the names that we were not able to resolve
    data.harm <- dplyr::filter(data.harm, !(is.na(scientificName) |is.na(db_taxon_higher_rank) | is.na(db_taxon_higher) ) )
    
    # join these data to the tax.dat data
    data <- dplyr::left_join(data, data.harm, by = c("row_id") )
    
    # add the special names back
    data <- dplyr::bind_rows(data, data.spec)
    
  } 
  
  # if we only have special names, then we only consider the special names and add additional columns for consistency
  else {
    
    data <- data.spec
    data[["scientificName"]] <- data[[clean.col]]
    data[["acceptedNameUsageID"]] <- NA
    data[["db_taxon_higher_rank"]] <- NA
    data[["db_taxon_higher"]] <- NA
    
  }
  
  # rename the life-stage column to life_stage
  data <- dplyr::rename(data, "life_stage" = dplyr::all_of(life_stage) )
  
  # remove the row_id column
  data <- dplyr::select(data, -row_id)
  
  # convert to a tibble
  data <- dplyr::as_tibble(data)
  
  return(data)
  
}


# generate some test data and run through the two functions
df.test <- 
  data.frame(taxon_name = c("Gammarus", "Daphnia", "Triops granitica", "Triops", 
                          "Simocephalus vetulus", "Turbellaria", "Nematoda"),
             Life_stage = c("adult", "adult", "adult", "adult",
                          "adult", "none", "none"),
             lat = rep(50.55, 7 ),
             lon = rep(4.98, 7),
             body_length_mm = rnorm(n = 7, mean = 10, sd = 1))
head(df.test)

# run the Clean_Taxon_Names() function
df.test1 <- Clean_Taxon_Names(data = df.test, 
                             target_taxon = "taxon_name", life_stage = "Life_stage",
                             database = "gbif")
head(df.test1)

# run the Get_Habitat_Data() function
df.test2 <- Get_Habitat_Data(data = df.test1, latitude_dd = "lat", longitude_dd = "lon")
View(df.test2)





