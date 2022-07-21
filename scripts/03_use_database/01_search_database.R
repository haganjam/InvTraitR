
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
#' @return data.frame with target taxon names and cleaned names for the chosen taxonomic backbone
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
    
    t1 <- (x[!is.na(x)] <= 90 & x[!is.na(x)] >= -90)
    
    t2 <- (y[!is.na(y)] <= 180 & y[!is.na(y)] >= -180)
    
    all( (t1 == TRUE) & (t2 == TRUE) )
    
  }
  
  assertthat::on_failure(test_6) <- function(call, env){
    
    paste0(deparse(call$x), " or ",deparse(call$x), " are too big/small to be valid decimal degrees")
    
  }
  
  assertthat::assert_that(test_6(x = data[[latitude_dd]], data[[longitude_dd]]))
  
  # create a data.frame with latitude and longitude columns
  lat.lon <- data.frame(longitude_dd = as.numeric(data[[longitude_dd]]),
                        latitude_dd = as.numeric(data[[latitude_dd]])) 
  head(lat.lon)
  
  # add a row id column
  lat.lon$row_id <- 1:nrow(lat.lon)
  
  # remove the NA values
  lat.lon2 <- lat.lon[!(is.na(lat.lon$longitude_dd) | is.na(lat.lon$latitude_dd )), ]
  
  # if there are no non-NA values then we don't get habitat data
  if( nrow(lat.lon2) == 0 ) {
    
    hab.dat <- dplyr::tibble(habitat_id = NA,
                             area_km2 = NA,
                             row_id = 1)
    
  } else {
    
    # convert this to a spatial points object
    sp.pts <- sp::SpatialPoints(lat.lon2[, c("longitude_dd", "latitude_dd")], proj4string = raster::crs(fw_map) )
    
    # check where pts overlap with freshwater ecoregion
    hab.dat <- sp::over(sp.pts, fw_map)
    
    # add row id to these data
    hab.dat$row_id <- lat.lon2$row_id
    names(hab.dat) <- c("habitat_id", "area_km2", "row_id")
    
    }
    
    # join to the original df to refill in the NAs
    hab.dat <- dplyr::full_join(lat.lon, hab.dat, by = "row_id")
    
    # arrange by row_id
    hab.dat <- dplyr::arrange(hab.dat, row_id)
    
    # remove the row_id column 
    hab.dat <- dplyr::select(hab.dat, -row_id)
    
    # join these data to the metadata
    hab.dat <- dplyr::left_join(hab.dat, fw_meta, by = "habitat_id")
    
    # bind the original data to the habitat data
    hab.dat <- dplyr::bind_cols(dplyr::select(data, -dplyr::all_of(c(latitude_dd, longitude_dd)) ), 
                                hab.dat)
    
    # convert to a tibble
    hab.dat <- dplyr::as_tibble(hab.dat)
    
    return(hab.dat)
  
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
  
  # make sure the correct packages are installed
  test_1 <- function(x) {
    
    all(c("bdc", "dplyr", "taxadb", "curl") %in% installed.packages()[,1])
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(c("bdc", "dplyr", "taxadb", "curl"), " must be installed for this function to work")
    
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
  
  # test if the target_taxon in the data object
  test_4 <- function(x, y) {
    
    assertthat::is.string(y) & ( y %in% names(x) )
    
  }
  
  assertthat::on_failure(test_4) <- function(call, env){
    
    paste0(deparse(call$y), " column is not a column in the supplied data object")
    
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
  
  # subset a name column from the data object
  name.dat <- data[target_taxon]
  
  # clean the names for typos etc. using the bdc_clean_names function
  clean.names <- bdc::bdc_clean_names(sci_names = name.dat[[target_taxon]], save_outputs = FALSE)
  
  # write some code to remove the output file
  unlink("Output", recursive=TRUE)
  
  # add the clean names to the data.frame
  clean.col <- paste("clean_", target_taxon, sep = "")
  name.dat[[clean.col]] <- clean.names$names_clean
  
  # check that all the non-missing names were cleaned
  test_6 <- function(x, y) {
    
    length(x[!is.na(x)]) == length(y[!is.na(y)])
    
  }
  
  assertthat::on_failure(test_6) <- function(call, env){
    
    "Length of clean names do not match length of original names"
    
  }
  
  assertthat::assert_that(test_6(x = name.dat[[target_taxon]], y = name.dat[[clean.col]]))
  
  # give each row an id
  name.dat$targ_no <- 1:nrow(name.dat)
  
  # update the database if there is a valid internet connection
  if (curl::has_internet()) {
    
    taxadb::td_create(provider = database,
                      overwrite = FALSE
    )
    
  }
  
  # subset out taxa with special names
  source(here::here("scripts/01_special_names_func.R"))
  spec.names <- special_taxon_names()
  name.dat.sp <- dplyr::filter(name.dat, (eval(parse(text = clean.col)) %in% spec.names) )
  
  # add a data.base column to the name.dat.sp object
  name.dat.sp[["db"]] <- "special"
  
  # remove the special names from the name.dat
  name.dat <- dplyr::filter(name.dat, !(eval(parse(text = clean.col)) %in% spec.names) )
  
  # harmonise the names to the gbif database
  harm.tax <- 
    bdc::bdc_query_names_taxadb(sci_name = name.dat[[clean.col]],
                                db = database,
                                rank_name = "Animalia",
                                rank = "kingdom"
    )
  
  # write some code to remove the output file
  unlink("Output", recursive=TRUE)
  
  # add columns for higher taxa, higher taxon ranks, db source and row_id
  
  # higher taxon rank
  harm.tax <- dplyr::mutate(harm.tax,
                            db_taxon_higher_rank = ifelse(is.na(order) & is.na(family), NA, 
                                                          ifelse(is.na(order) & !is.na(family), "family", "order") ) )
  # higher taxon name
  harm.tax <- dplyr::mutate(harm.tax,
                            db_taxon_higher = ifelse(is.na(order) & is.na(family), NA, 
                                                     ifelse(is.na(order), family, order) ) )                          
  
  # row_id
  harm.tax$targ_no <- 1:nrow(harm.tax)
  
  # select the relevant columns
  harm.tax <- harm.tax[,c("targ_no", "scientificName", "acceptedNameUsageID", "db_taxon_higher_rank", "db_taxon_higher")]
  
  # remove the names that we were not able to resolve
  harm.tax <- dplyr::filter(harm.tax, !(is.na(scientificName) |is.na(db_taxon_higher_rank) | is.na(db_taxon_higher) ) )
  
  # add the database to the name.dat
  name.dat$db <- database
  
  # join these data to the tax.dat data
  name.dat <- dplyr::left_join(name.dat, harm.tax, by = c("targ_no") )
  
  # add the special names back
  name.dat <- dplyr::bind_rows(name.dat, name.dat.sp)
  
  # make sure the targ_no column is a character
  name.dat <- dplyr::mutate(name.dat, targ_no = as.character(targ_no))
  
  # rename the life-stage column to life_stage
  data <- dplyr::rename(data, "life_stage" = dplyr::all_of(life_stage) )
  
  # join the rest of the input data
  name.dat <- dplyr::full_join(name.dat, data, by = c(target_taxon, "targ_no") )
  
  # remove the targ_no column
  name.dat <- dplyr::select(name.dat, -targ_no)
  
  # convert to a tibble
  name.dat <- dplyr::as_tibble(name.dat)
  
  return(name.dat)
  
}

#'
#' @title Select_Traits()
#' 
#' @description Get taxonomic distances of target names relative to the taxa databases
#' 
#' @details This function searches the relevant trait or equation database for the best
#' matching trait or equation for a given target name based on three criteria: taxonomic distance
#' life-stage match and habitat match
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param input - input data.frame exported from Get_Habitat_Data() and Clean_Taxon_Names() function
#' @param target_taxon - character string with the column name containing the taxon names
#' @param max_tax_dist - maximum taxonomic distance acceptable between the target and the taxa in the database (default = 3)
#' @param trait - trait to be searched for (default = "equation")
#' @param workflow - "workflow1" or "workflow2" (default)
#' @param gen_sp_dist - taxonomic distance between a genus and a species (default = 0.5)
#' 
#' @return data.frame with target taxon names and cleaned names for the chosen taxonomic backbone
#' 

Select_Traits <- function(input, target_taxon, 
                          max_tax_dist = 3, trait = "equation", 
                          workflow = "workflow2", 
                          gen_sp_dist = 0.5) {
  
  # make sure the correct packages are installed
  test_1 <- function(x) {
    
    all(c("igraph", "dplyr") %in% installed.packages()[,1])
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(c("igraph", "dplyr"), " must be installed for this function to work")
    
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
  if (!exists("fw")) {
    fw <- readRDS(file = here::here("database/freshwater_ecoregion_data.rds"))
  }
  
  # split the input data.frame into a list
  input.list <- split(input, 1:nrow(input))
  
  # for each entry in the input.list, select appropriate traits
  
  output <- 
    
    lapply(input.list, function(data) {
      
      
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
          
          return( dplyr::tibble(targ.scientificName = data[[paste0("clean_", target_taxon)]]) )
          
        }
        
        # if the higher taxon is not in the database then return an NA
        if (all((names(htm_db) == data[["db_taxon_higher"]]) == FALSE)) {
          
          return( dplyr::tibble(targ.scientificName = data[[paste0("clean_", target_taxon)]]) )
          
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
      
    })
  
    output <- dplyr::bind_rows(output, .id = "targ_no")
    
    # if workflow1, then output all equations within the correct taxonomic distance
    if (workflow == "workflow1") {
      
      x <- dplyr::rename(cleaned.names, targ.scientificName = paste0("clean_", target_taxon) )
      y <- dplyr::full_join(x, output, by = "targ.scientificName")
      names(y)[names(y) == "id"] <- paste0(trait, "_id")
      output <- dplyr::left_join(y, trait_db, by = paste0(trait, "_id"))
      
    }
    
    return(output)
  
}

#'
#' @title Get_Trait_From_Taxon()
#' 
#' @description Get taxonomic distances of target names relative to the taxa databases
#' 
#' @details This function takes a data.frame with target_taxa and uses the bdc package
#' to clean and harmonise the names.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param input_data - data.frame with a columns: target taxon name, life-stage, body size, 
#' and latitude and longitude (decimal degrees)
#' @param target_taxon - character string with the column name containing the taxon names
#' @param life_stage - character string with the column name containing the life_stages
#' @param body_size - character string with the column name containing the body size
#' @param latitude_dd - character string with the column name containing the latitude (decimal degrees)
#' @param longitude_dd - character string with the column name containing the longitude (decimal degrees)
#' @param trait - name of the trait to search (only "equation" is currently supported)
#' @param workflow - "workflow1" or "workflow2" (default)
#' @param max_tax_dist - maximum suitable taxonomic distance
#' @param gen_sp_dist - taxonomic distance between a genus and a species
#' 
#' @return data.frame with outputted trait information
#' 

Get_Trait_From_Taxon <- function(input_data,
                                 target_taxon, life_stage, body_size,
                                 latitude_dd, longitude_dd,
                                 trait = "equation",
                                 workflow = "workflow2", 
                                 max_tax_dist = 3, gen_sp_dist = 0.5
                                 ) {
  
  # make sure the correct packages are installed
  test_1 <- function(x) {
    
    all(c("dplyr") %in% installed.packages()[,1])
    
  }
  
  assertthat::on_failure(test_1) <- function(call, env){
    
    paste0(c("dplyr"), " must be installed for this function to work")
    
  }
  
  assertthat::assert_that(test_1())
  
  # make sure the max_tax_dist argument is a number
  test_2 <- function(x, y) {
    
    all(x %in% names(y))
    
  }
  
  assertthat::on_failure(test_2) <- function(call, env){
    
    paste0(deparse(call$x), " are not all present in the input_data object")
    
  }
  
  assertthat::assert_that(test_2( x = c(target_taxon, life_stage, body_size, latitude_dd, longitude_dd),
                                  y = input_data))
  
  
  # get unique data points
  unique_data <- dplyr::distinct( dplyr::select(input_data, all_of(c(target_taxon, life_stage, latitude_dd, longitude_dd)) ) )
  
  # get the habitat data
  hab.dat <- Get_Habitat_Data(data = unique_data, latitude_dd = latitude_dd, longitude_dd = longitude_dd)
  hab.dat[["targ_no"]] <- as.character(1:nrow(hab.dat))
  
  # make a vector of taxon databases
  tax.vector <- c("gbif", "itis", "col")
  
  trait.dat <- 
    
    lapply(tax.vector, function(tax_database) {
      
      # clean the names and place in one of the three taxonomic backbones
      cleaned.names <- Clean_Taxon_Names(data = hab.dat, 
                                         target_taxon = target_taxon,
                                         life_stage = life_stage, database = tax_database)
      
      # select the traits based on taxonomic distance, life-stage and habitat
      selected.traits <- Select_Traits(input = cleaned.names,
                                       target_taxon = target_taxon,
                                       max_tax_dist = max_tax_dist, 
                                       trait = trait, 
                                       workflow = workflow,
                                       gen_sp_dist = gen_sp_dist)
      
      return(selected.traits)
      
    })
  
  names(trait.dat) <- tax.vector
  trait.dat <- dplyr::bind_rows(trait.dat,  .id = "tax_database")
  
  if (workflow == "workflow1") {
    
    return(trait.dat)
    
  }
  
  # how to select traits from the list
  trait.dat <- 
    
    lapply(split(trait.dat, trait.dat$targ_no), function(x){
      
      # if all the taxonomic distances are NA then return an NA
      if( all(is.na(x$tax_distance)) ) {
        
        return(x[sample(1:nrow(x), 1), c("targ_no") ])
        
      } else {
        
        # remove NAs from tax_distance
        y <- x[!is.na(x$tax_distance), ]
        
        # select the minimum taxonomic distance
        z <- y[ dplyr::near(y$tax_distance, min(y$tax_distance, na.rm = TRUE)), ]
        
        # calculate habitat match score
        u <- apply(z[,c("realm_match", "maj_hab_match", "ecoregion_match")], 1, sum)
        
        # subset the entries with the highest habitat match score
        if (all(!is.na(u))) {
          
          z <- z[u == max(u), ]
          
        }
        
        # if there are still multiple equations, choose randomly
        z <- z[sample(1:nrow(z), 1), ] 
        
        return(z)
        
      }
      
    })
  
  output <- dplyr::full_join(hab.dat, dplyr::bind_rows(trait.dat), by = "targ_no")
  output <- dplyr::select(output, -targ_no)
 
  # add other information to the input data after we have got unique data points
  output <- 
    
    dplyr::left_join(dplyr::rename(input_data, 
                                   latitude_dd = dplyr::all_of(latitude_dd), 
                                   longitude_dd = dplyr::all_of(longitude_dd)), 
                     output, 
                     by = c(target_taxon, life_stage, "latitude_dd", "longitude_dd")
    )
  
  # evaluate the trait i.e. get the trait value
  if (!exists(paste0(trait, "_db"))) {
    
    assign( paste0(trait, "_db"),
            readRDS(file = paste0(here::here("database"), "/", trait, "_database.rds")))
    
  }
  
  # assign the object to trait_db
  trait_db <- get(paste0(trait, "_db"))
  
  # assign the relevant trait values
  if (trait == "equation") {
    
    weight_mg <- vector(length = nrow(output))
    default_bs <- vector(length = nrow(output))
    flags <- vector(length = nrow(output))
    
    for(i in 1:nrow(output)) {
      
      id <-  output[i,][["id"]]
      trait_sel <- trait_db[trait_db[["equation_id"]] == id, ]
      var1 <- output[i, ][[body_size]]
      
      # if the body size variable is missing (i.e. NA), replace with midpoint
      if (is.na(var1)) {
        
        var1 <- sum(trait_sel[["body_size_min"]], trait_sel[["body_size_max"]])/2
        default_bs[i] <- var1
        
        flags[i] <- "used midpoint from equation"
        weight_mg[i] <- eval(parse(text = trait_sel[["equation"]]))
        
      } else {
        
        default_bs[i] <- NA
        
        range_min <- (var1 > trait_sel[["body_size_min"]])
        range_max <- (var1 < trait_sel[["body_size_max"]])
        
        flags[i] <- ifelse(any(c(range_min, range_max) == FALSE), 
                           "inputted length data is beyond the range used to fit the equation",
                           NA)
        
        weight_mg[i] <- eval(parse(text = trait_sel[["equation"]]))
        
      }
      
    }
    
    output[["default_body_size"]] <- default_bs
    output[["weight_mg"]] <- weight_mg
    output[["flags"]] <- flags
    
  } else {
    
    output <- dplyr::left_join(output, trait_db[, c(paste0(trait, "_id"), trait)])
    
  }
  
  # convert the output to a tibble
  output <- dplyr::as_tibble(output)
  
  return(output)
  
}

### END
