
#'
#' @title Extract_Genus()
#' 
#' @description If a binomial taxa is supplied, extract the genus name only
#' 
#' @details This method works almost exclusively with genera and not species. This is because
#' once the genus is known, the relationship among the species within that genus is known. It
#' makes the method much more computationally efficient. This function is used to then extract
#' the genus from a species name. The genus name is used to place the species in the correct
#' taxonomic framework 
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param binomial - binomial character string separated by a space (e.g. "Loxodonta loxodonta")
#' 
#' @return string with the genus name
#' 

Extract_Genus <- function(binomial) {
  
  # test if the input is a string
  test_1 <- function(x) {
    
    assertthat::is.string(x)
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(deparse(call$x), " is not a character string")
    
  }
  
  assertthat::assert_that(test_1(x = binomial))
  
  # test if the input is of length 1
  test_2 <- function(x) {
    
    assertthat::are_equal(length(x), 1)
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " is not of length = 1")
    
  }
  
  assertthat::assert_that(test_2(x = binomial))
  
  # check if input has any special characters
  test_3 <- function(x) {
    
    !grepl(pattern = '[[:punct:]]', x = x)
    
  }
  
  assertthat::on_failure(test_3) <- function(call, env){
    
    paste0(deparse(call$x), " contains special characters")
    
  }
  
  assertthat::assert_that(test_3(x = binomial))
  
  # split the binomial into separate parts
  binomial.1st <- unlist( strsplit(x = binomial, split = " ", fixed = TRUE) )
  
  # calculate the length of the split object
  binomial.l <- length(binomial.1st)
  
  # if the resulting object has a length of greater than 1, then extract first element
  if (binomial.l > 1) {
    
    binomial <- binomial.1st[1]
    
  } 
  
  # add a word count attribute
  attr( binomial, "n") <- binomial.l
  
  # return the modified name
  return(binomial)
  
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
    data[["scientificName"]] <- NA
    data[["acceptedNameUsageID"]] <- NA
    data[["db_taxon_higher_rank"]] <- NA
    data[["db_taxon_higher"]] <- NA
    
  }
  
  # remove the row_id column
  data <- dplyr::select(data, -row_id)
  
  # convert to a tibble
  data <- dplyr::as_tibble(data)
  
  return(data)
  
}

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
  
  # make a vector of packages
  pack_vec <- c("dplyr", "sp", "raster")
  
  # make sure the correct packages are installed
  test_1 <- function(x) {
    
    all( pack_vec %in% installed.packages()[,1])
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(pack_vec, " must be installed for this function to work")
    
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

# generate some test data and run through the two functions
df.test <- 
  data.frame(taxon_name = c("Gammarus", "Daphnia", "Triops granitica", "Triops", 
                          "Simocephalus vetulus", "Turbellaria", "Nematoda"),
             Life_stage = c("adult", "adult", "adult", "adult",
                          "adult", "none", "none"),
             lat = rep(50.55, 7 ),
             lon = rep(4.98, 7))
head(df.test)

# run the Clean_Taxon_Names() function
df.test1 <- Clean_Taxon_Names(data = df.test, 
                              target_taxon = "taxon_name", life_stage = "Life_stage",
                              database = "gbif")
head(df.test1)

# run the Get_Habitat_Data() function
df.test2 <- Get_Habitat_Data(data = df.test1, latitude_dd = "lat", longitude_dd = "lon")
View(df.test2)

data = df.test2
target_taxon <- "target_name"
life_stage <- "Life_stage"
max_tax_dist = 3
trait = "equation"
workflow = "workflow2" 
gen_sp_dist = 0.5

# make a vector of packages
pack_vec <- c("igraph", "dplyr")

# make sure the correct packages are installed
test_1 <- function(x) {
  
  all(pack_vec %in% installed.packages()[,1])
  
}

assertthat::on_failure(test_1) <- function(call, env){
  
  paste0(pack_vec, " must be installed for this function to work")
  
}

assertthat::assert_that(test_1())

# make sure the max_tax_dist argument is a number
test_2 <- function(x) {
  
  assertthat::is.number(x) & (x >= 0)
  
}

assertthat::on_failure(test_2) <- function(call, env){
  
  paste0(deparse(call$x), " must be a number that is either greater than or equal to zero")
  
}

assertthat::assert_that(test_2(max_tax_dist))
assertthat::assert_that(test_2(gen_sp_dist))

# make sure the max_tax_dist argument is a number
test_3 <- function(x) {
  
  trait %in% c("equation", paste0("trait", 1:10))
  
}

assertthat::on_failure(test_3) <- function(call, env){
  
  paste0(deparse(call$x), " is not a valid trait or equation, see documentation")
  
}

assertthat::assert_that(test_3(trait))

# load the igraph package so that igraph objects can be manipulated
library(igraph)

# evaluate the trait i.e. get the trait value
if (!exists(paste0(trait, "_db"))) {
  
  assign( paste0(trait, "_db"),
          readRDS(file = paste0(here::here("database"), "/", trait, "_database.rds")))
  
}

# assign the object to trait_db
trait_db <- get(paste0(trait, "_db"))

# load the habitat database
if (!exists("hab_db")) {
  hab_db <- readRDS(file = here::here("database/freshwater_ecoregion_data.rds"))
}

# add an identifier for the individual inputs
data[["row_id"]] <- 1:nrow(data)

# split the input data.frame into a list
data.list <- split(data, data[["row_id"]])
data.list

# for each entry in the input.list, select appropriate traits

# output <- 
  
  # lapply(data.list, function(input) {

# there are essentially thee options:
# 1. scientificName is not NA which means that name is recognised by the taxonomic backbone
# 2. scientificName is NA and database is either gbif, itis or col which means that the name is not recognised by the taxonomic backbone
# 3. scientificName is NA but the database is special which means it is a special name

# if scientificName is not NA

input <- data.list[[1]]

if ( !is.na( input[["scientificName"]] ) & !is.na( input[["db_taxon_higher"]] ) ) {TRUE}



# load the taxon matrices
if (!exists("htm_db")) {
  htm_db <- readRDS(file = paste0(here::here("database"), "/", input[["db"]], "_higher_taxon_matrices.rds"))
}

# load the taxon database
if (!exists("td_db")) {
  td_db <- readRDS(file = paste0(here::here("database"), "/", input[["db"]], "_taxon_database.rds"))
}  

# taxon matrix
htm <- htm_db[ ( names(htm_db) == input[["db_taxon_higher"]] ) ][[1]]

# extract vertices
v.x <- V(htm)

# extract the equation entries from the taxon database
td <- td_db[td_db$database == trait,]

# extract the entries from the equation taxon database matching the target higher taxon
td <- td[td$db_taxon_higher == input[["db_taxon_higher"]], ]

# extract the target.name
target.name <- input[["scientificName"]]

# taxonomic distance
dist.df <- 
  
  mapply( function(db.name, id) {

    # extract genus for species-level names
    target.name2 <- Extract_Genus(target.name)
    db.name2 <- Extract_Genus(db.name)
    
    tax.dist <- 
      
      igraph::distances(htm, 
                        v.x[which(attr(v.x, "names") == target.name2)],
                        v.x[which(attr(v.x, "names") == db.name2)],
                        mode = c("all"),
                        algorithm = c("bellman-ford")
                        
      )
    
    # if length is zero then the distance is zero
    if(length(tax.dist) == 0) {
      
      tax.dist <- NA
      
    } else {
      
      tax.dist <- tax.dist[[1]]
      
    }
    
    # extra distance for species level: gen_sp_dist argument
    sp.l <- sum( ifelse( c(attr(target.name2, "n"), attr(db.name2, "n")) > 1,  gen_sp_dist, 0))
    
    dist.df <- 
      dplyr::tibble(scientificName = target.name,
                    db.scientificName = db.name,
                    trait_out = trait,
                    id = id,
                    tax_distance = tax.dist + sp.l
      )
    
    return(dist.df)
    
  }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE )

# bind into a data.frame
dist.df <- dplyr::bind_rows(dist.df)


#### I WAS HERE



# if the taxonomic distances are too big then we output NA
if (min(tax.dist.df[["tax_distance"]]) > max_tax_dist ) {
  
  return(dplyr::tibble(targ.scientificName = data[["scientificName"]]))
  
}

# remove the rows where tax_distance is too large
tax.dist.df <- dplyr::filter(tax.dist.df, tax_distance <= max_tax_dist)

}

# if workflow1, then output all equations within the correct taxonomic distance
if (workflow == "workflow1") {
  
  return(tax.dist.df)
  
}









# else if ( is.na(input[["scientificName]]) & (input[["db]] == "special) ) {}

# else if ( is.na(input[["scientificName"]] & ( any( input[["db"]] == c("gbif", "itis", "col") ) )) ) {}
    


    if (data[["db"]] == "special") {
      
      trait_spec <- dplyr::filter(trait_db, db_taxon == data[[paste0("clean_", target_taxon)]] )
      
      tax.dist.df <- 
        dplyr::tibble(targ.scientificName = data[[paste0("clean_", target_taxon)]],
                      db.scientificName = if(nrow(trait_spec) == 0) { NA } else { trait_spec[["db_taxon"]] },
                      trait_out = trait,
                      id = if(nrow(trait_spec) == 0) { NA } else { trait_spec[[paste0(trait, "_id")]] },
                      tax_distance = NA)
      
    } else {
      
      # load the taxon matrices
      if (!exists("htm_db")) {
        htm_db <- readRDS(file = paste0(here::here("database"), "/", data[["db"]], "_higher_taxon_matrices.rds"))
      }
      
      # load the taxon database
      if (!exists("td_db")) {
        td_db <- readRDS(file = paste0(here::here("database"), "/", data[["db"]], "_taxon_database.rds"))
      }  
      
      # if the input is NA then we return an NA
      if (any( is.na(c(data[["scientificName"]], data[["db_taxon_higher"]])) ) ) {
        
        return( dplyr::tibble(targ.scientificName = data[["scientificName"]]) )
        
      }
      
      # if the higher taxon is not in the database then return an NA
      if (all((names(htm_db) == data[["db_taxon_higher"]]) == FALSE)) {
        
        return( dplyr::tibble(targ.scientificName = data[["scientificName"]]) )
        
      }
      
      # taxon matrix
      htm <- htm_db[(names(htm_db) == data[["db_taxon_higher"]])][[1]]
      
      # extract vertices
      v.x <- V(htm)
      
      # extract the taxon database
      td <- td_db[td_db$database == trait,]
      td <- td[td$db_taxon_higher == data[["db_taxon_higher"]], ]
      
      # extract the target.name
      targ.name <- data[["scientificName"]]
      
      # taxonomic distance
      
      tax.dist.df <- 
        mapply( function(db.name, id) {
          
          # extract genus for species-level names
          targ.name2 <- Extract_Genus(targ.name)
          db.name2 <- Extract_Genus(db.name)
          
          tax.dist <- 
            igraph::distances(htm, 
                              v.x[which(attr(v.x, "names") == targ.name2)],
                              v.x[which(attr(v.x, "names") == db.name2)],
                              mode = c("all"),
                              algorithm = c("bellman-ford")
            )
          
          # if length is zero then the distance is zero
          if(length(tax.dist) == 0) {
            
            tax.dist <- 0
            
          } else {
            
            tax.dist <- tax.dist[[1]]
            
          }
          
          # extra distance for species level: gen_sp_dist argument
          sp.l <- sum( ifelse(c(attr(targ.name2, "n"), attr(db.name2, "n")) > 1,  gen_sp_dist, 0))
          
          dist.df <- 
            dplyr::tibble(targ.scientificName = targ.name,
                          db.scientificName = db.name,
                          trait_out = trait,
                          id = id,
                          tax_distance = tax.dist + sp.l
            )
          
          return(dist.df)
          
        }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE)
      
      # bind into a data.frame
      tax.dist.df <- dplyr::bind_rows(tax.dist.df)
      
      # if the taxonomic distances are too big then we output NA
      if (min(tax.dist.df[["tax_distance"]]) > max_tax_dist ) {
        
        return(dplyr::tibble(targ.scientificName = data[["scientificName"]]))
        
      }
      
      # remove the rows where tax_distance is too large
      tax.dist.df <- dplyr::filter(tax.dist.df, tax_distance <= max_tax_dist)
      
    }
    
    # if workflow1, then output all equations within the correct taxonomic distance
    if (workflow == "workflow1") {
      
      return(tax.dist.df)
      
    }
    
    # life-stage data
    
    # subset the equations that have the correct taxonomic distance
    trait_sel <- trait_db[trait_db[[paste0(trait, "_id")]] %in% tax.dist.df$id, ]
    
    # subset the equations with the correct life-stage
    trait_sel <- dplyr::filter(trait_sel, db_life_stage == data[["life_stage"]])
    
    # subset the max_tax_dist
    tax.dist.df <- dplyr::filter(tax.dist.df, id %in% trait_sel[[paste0(trait, "_id")]])
    tax.dist.df[["life_stage_match"]] <- TRUE
    
    # from the equations within the relevant maximum taxonomic distance and with the correct life-stage
    # choose the equation with the lowest taxonomic distance
    if(any(!is.na(tax.dist.df$tax_distance)) ) {
      
      tax.dist.df <- tax.dist.df[dplyr::near(tax.dist.df$tax_distance, min(tax.dist.df$tax_distance)), ]
      
    }
    
    # habitat data
    
    # get the relevant trait
    fw_sel <- dplyr::filter(fw, database == trait)
    
    # get the equations with the correct taxonomic distance and life-stage
    fw_sel <- dplyr::filter(fw_sel, id %in% tax.dist.df[["id"]])
    
    # filter progressively based on habitat
    tax.dist.df[["realm_match"]] <- (fw_sel[["realm"]] == data[["realm"]])
    tax.dist.df[["maj_hab_match"]] <- (fw_sel[["major_habitat_type"]] == data[["major_habitat_type"]])
    tax.dist.df[["ecoregion_match"]] <- (fw_sel[["ecoregion"]] == data[["ecoregion"]])
    
    # add them up
    hab.sim <- tax.dist.df[["realm_match"]] + tax.dist.df[["maj_hab_match"]] + tax.dist.df[["ecoregion_match"]]
    
    # select the best matching habitat types
    if( any(!is.na(hab.sim)) ) {
      
      # select the best matching habitat
      tax.dist.df <- tax.dist.df[dplyr::near(hab.sim, max(hab.sim)), ]
      tax.dist.df[["habitat_flag"]] <- "none"
      
    } else {
      
      # give a warning if none of the habitats match
      tax.dist.df[["habitat_flag"]] <- as.character(ifelse(hab.sim == 0, "different realm, major habitat type and ecoregion", "none"))
      
    }
    
    # make sure the db.scientificName column is a character variable
    tax.dist.df[["db.scientificName"]] <- as.character(tax.dist.df[["db.scientificName"]])
    
    return(tax.dist.df)
    
  # })

output <- dplyr::bind_rows(output, .id = "targ_no")

# join the input data
in_out_join <- dplyr::full_join(input, output, by = "targ_no")

# reorder the columns
in_out_join <- 
  in_out_join %>%
  dplyr::select(targ_no, dplyr::all_of(names(in_out_join)[names(in_out_join) != "targ_no"] ))

# if workflow1, then output all equations within the correct taxonomic distance
if (workflow == "workflow1") {
  
  names(in_out_join)[names(in_out_join) == "id"] <- paste0(trait, "_id")
  
  # join the equations
  in_out_join <- dplyr::left_join(in_out_join, trait_db, by = paste0(trait, "_id"))
  
  return(in_out_join)
  
} else if (workflow == "workflow2") {
  
  return(in_out_join)
  
} else {
  
  stop("Choose an appropriate workflow: workflow1 or workflow2")
  
}



