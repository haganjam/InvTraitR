
# load relevant libraries
library(here)

# check for the correct packages:
source(here("scripts/create_database/01_version_package_warnings.R"))

# load relevant functions
source(here("scripts/functions/01_get_taxon_id_function.R"))
source(here("scripts/create_database/03_gbif_downstream_function.R"))
source(here("scripts/create_database/04_itis_downstream_function.R"))


#' @title get_taxon_order()
#' 
#' @description Function to take any taxon name and output the order name
#' 
#' @details For a given entry in the equation or length database, the function
#' will take the taxon name in the database and find the order. It then wraps this
#' order information with other relevant information like the equation or length ID
#' and the life.stage. Currently, only "itis" is supported but there are plans to expand
#' this to the gbif database as well.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param equ.name.input - name in the equation or length database
#' @param equ.id - equation or length id in the databases
#' @param data.base - "itis" and "gbif" are supported
#' @param life.stage - taxon life stage in the equation or length database
#' 
#' @return list with two elements:
#' 1. taxlist - list with the taxon name, whether it is a suitable equation, the database
#' the rank.name, the order name and the life-stage
#' 2. synonymns - list with a vector of synonymns for that taxon name
#' 

get_taxon_order <- function(name.input, data.base) {
  
  # create an output list
  taxlist <- list(name = name.input,
                  accepted_name = NA,
                  name_rank = NA,
                  suitable = FALSE,
                  database = data.base,
                  higher_rank = NA,
                  higher_name = NA)
  
  # if the input name is a species then extract the genus
  name.in <- extract_genus(binomial = name.input)
  
  # get the taxon_id from the correct database
  taxon_id <- get_taxon_id(database_function = data.base, taxon_name = name.in, ask_or_not = FALSE, tries = 5)
  
  # if the taxon_id is not found in the database, then stop the function
  
  # check if is of length == 0
  if(length(taxon_id) == 0) {
    
    taxon_id <- NA
    
  }
  
  # then, check if it is NA (both are possible)
  if (is.na(taxon_id)) {
    
    return(taxlist)
    stop("Could not find this taxon name in the database")
    
  }
  
  # get the upwards classification for the taxon_id
  y.c <- classification(taxon_id)[[1]]
  y.c$ranknum <- 1:nrow(y.c)
  
  # check that the rank is Animalia otherwise return an error
  if (y.c[y.c$rank == "kingdom",]$name != "Animalia" ) {
    
    return(taxlist)
    stop("Classification is incorrect as all taxa are in the Animalia kingdom")
    
  }
  
  # add the accepted name to the taxlist database
  taxlist$accepted_name <- y.c[nrow(y.c),][["name"]] 
  
  # get the taxonomic rank name for the taxon_id
  if (length(unlist( strsplit(x = name.input, split = " ", fixed = TRUE) )) > 1) {
    
    taxlist$name_rank <- "species"
    
  } else { 
    
    taxlist$name_rank <- y.c[nrow(y.c),][["rank"]] 
    
  }
  
  # check if the classification has an order
  above_order <- c("kingdom", "subkingdom", "infrakingdom", 
                   "superphylum", "phylum", "subphylum", "infraphylum", 
                   "superclass", "class", "subclass", "infraclass", "superorder"
  )
  
  # if the taxon_id is above the rank of order, then the equation is unsuitable
  if( taxlist$name_rank %in% above_order ) {
    
    taxlist$suitable <- FALSE
    
  } else {
    
    taxlist$suitable <- TRUE
    
  }
  
  # if order is in the higher classification then we add that
  if( "order" %in% y.c$rank) {
    
    # specify order as the higher rank and add the name
    taxlist$higher_rank <- "order"
    
    # get the order name and rank number for the taxon_id
    taxlist$higher_name <- y.c[y.c$rank == "order", ][["name"]]
    
  } else {
    
    # get the highest rank leftover that isn't above order
    max_rank <- y.c$rank[!(y.c$rank %in% above_order)][1]
    
    # add that as higher rank
    taxlist$higher_rank <- max_rank
    
    # add the higher rank name
    taxlist$higher_name <- y.c[y.c$rank == max_rank, ][["name"]]
    
  }
  
  return(taxlist)
  
}

#'
#' @title get_taxon_distance()
#' 
#' @description Get a taxonomic distance matrix of all taxa from a given order
#' 
#' @details This function will take an order name and get all the taxa that are downstream
#' from that order and process it into a taxon matrix.
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param higher.tax.name - name of the higher taxa to get downstream taxa from
#' @param higher.tax.rank - rank of the higher.tax.name
#' @param data.base - name of the taxonomic database (only "itis" is currently supported)
#' 
#' @return igraph object with the taxonomic distance matrix
#' 

get_taxon_distance <- function(higher.tax.name, higher.tax.rank = NA, data.base = "itis") {
  
  # get taxon id: database specific function i.e. get the taxon ID of the order of the taxa from the database
  ord.id <- get_taxon_id(database_function = data.base, 
                         taxon_name =  higher.tax.name, 
                         ask_or_not = FALSE, tries = 5)
  ord.id <- ord.id[[1]]
  
  # create an error dmat output
  dmat <- c(list(higher.tax = higher.tax.name),
            list(tax_available = FALSE),
            list(tax_distance = NA,
                 tax_names = NA))
  
  # if the taxon_id is not found in the database, then stop the function
  if (is.na(ord.id)) {
    
    return(dmat)
    stop("Could not find this order name in the database")
    
  }
  
  # get a distance matrix and list of taxa downstream of the order
  if (data.base == "gbif") {
    
    dmat <- downstream_gbif(ord.id = ord.id, 
                            ord.name = higher.tax.name,
                            higher.tax.rank = higher.tax.rank)
    
  } else if (data.base == "itis") {
    
    dmat <- downstream_itis(ord.id = ord.id, ord.name = higher.tax.name)
    
  } else {
    
    stop("Specify an appropriate database: itis or gbif")
    
  }
  
  if (identical(dmat, NULL)) {
    
    return(dmat)
    stop("Could not find this order name in the database")
    
  }
  
  dmat <- c(list(higher.tax = higher.tax.name),
            list(tax_available = TRUE),
            dmat)
  
  return(dmat)
  
}

### END
