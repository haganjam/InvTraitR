
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

#'
#' @title Select_Traits_Tax_Dist()
#' 
#' @description Get taxonomic distances of target names relative to the taxa databases
#' 
#' @details This function searches the relevant trait or equation database for the best
#' matching trait or equation for a given target name based on three criteria: taxonomic distance
#' life-stage match and habitat match
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param data - input data.frame exported from Get_Habitat_Data() and Clean_Taxon_Names() function
#' @param target_taxon - character string with the column name containing the taxon names
#' @param max_tax_dist - maximum taxonomic distance acceptable between the target and the taxa in the database (default = 3)
#' @param trait - trait to be searched for (default = "equation")
#' @param gen_sp_dist - taxonomic distance between a genus and a species (default = 0.5)
#' 
#' @return tibble of the input data with traits or equations within the maximum taxonomic distance
#' 

Select_Traits_Tax_Dist <- function(data, 
                                   target_taxon, 
                                   max_tax_dist = 3, trait = "equation",
                                   gen_sp_dist = 0.5) {
  
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
  
  # load the trait data
  if (!exists(paste0(trait, "_db"))) {
    
    assign( paste0(trait, "_db"),
            readRDS(file = paste0(here::here("database"), "/", trait, "_database.rds")))
    
  }
  
  # assign the object to trait_db
  trait_db <- get(paste0(trait, "_db"))
  
  # load the taxon matrices
  if (!exists("htm_db")) {
    htm_db <- readRDS(file = paste0(here::here("database"), "/", input[["db"]], "_higher_taxon_matrices.rds"))
  }
  
  # load the taxon database
  if (!exists("td_db")) {
    td_db <- readRDS(file = paste0(here::here("database"), "/", input[["db"]], "_taxon_database.rds"))
  } 
  
  # add an identifier for the individual inputs
  data[["row_id"]] <- 1:nrow(data)
  
  # split the input data.frame into a list
  data.list <- split(data, data[["row_id"]])
  
  # for each entry in the input.list, select appropriate traits
  output <- 
    
    lapply(data.list, function(input) {
      
      if ( !is.na( input[["scientificName"]] ) & !is.na( input[["db_taxon_higher"]] ) & ( input[["db_taxon_higher"]] %in% names(htm_db) ) ) {
        
        # taxon matrix
        htm <- htm_db[ ( names(htm_db) == input[["db_taxon_higher"]] ) ][[1]]
        
        # extract vertices
        v.x <- V(htm)
        
        # extract the equation entries from the taxon database
        td <- td_db[td_db$database == trait, ]
        
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
              dplyr::tibble(db.scientificName = db.name,
                            trait_out = trait,
                            id = id,
                            tax_distance = tax.dist + sp.l
              )
            
            return(dist.df)
            
          }, td[["scientificName"]], td[["id"]], SIMPLIFY = FALSE )
        
        # bind into a data.frame
        dist.df <- dplyr::bind_rows(dist.df)
        
        # remove the rows where tax_distance is too large
        dist.df <- dplyr::filter(dist.df, tax_distance <= max_tax_dist)
        
      } else if ( is.na(input[["scientificName"]]) & (input[["db"]] == "special") ) {
        
        # get row_id's from trait database matching the special names
        row_id <- which(trait_db[["db_taxon"]] == input[[paste0("clean_", target_taxon)]])
        
        # check if there are rows that outputted and if not return NA
        x <- if(length(row_id) == 0) {NA} else {trait_db[row_id, ][["db_taxon"]]}
        y <- if(length(row_id) == 0) {NA} else {trait_db[row_id, ][[paste0(trait, "_id")]]}
        
        # pull this into a data.frame
        dist.df <- dplyr::tibble(db.scientificName = x,
                                 trait_out = trait,
                                 id = y,
                                 tax_distance = NA)
        
      } else {
        
        message("Input data do not match any of the conditions for selecting appropriate traits/equations")
        message(paste0("Returning NA for ", input[[paste0("clean_", target_taxon)]]) )
        
        dist.df <- dplyr::tibble(db.scientificName = NA,
                                 trait_out = trait,
                                 id = NA,
                                 tax_distance = NA
        )
        
      }
      
      # add metadata
      dist.df <- dplyr::bind_cols(input, dist.df)
      
      return(dist.df)
      
    } ) 
  
  return(output)
  
  }


trait_db
View(z1[[1]])

# load the habitat database
if (!exists("hab_db")) {
  hab_db <- readRDS(file = here::here("database/freshwater_ecoregion_data.rds"))
}

# load the trait data
if (!exists(paste0(trait, "_db"))) {
  
  assign( paste0(trait, "_db"),
          readRDS(file = paste0(here::here("database"), "/", trait, "_database.rds")))
  
}

# assign the object to trait_db
trait_db <- get(paste0(trait, "_db"))


# add database life-stage information
db_life_stage <- trait_db[trait_db[[paste0(trait, "_id")]] %in% dist.df[["id"]], ][["db_life_stage"]]
dist.df[["db_life_stage"]] <- if(length(db_life_stage) == 0) {NA} else {db_life_stage}

View(z1[[1]])


