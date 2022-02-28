
# Get taxonomic information for each equation in the database

# load relevant libraries
library(here)
library(igraph)
library(taxize)
library(dplyr)

# check for the correct packages:
source(here("scripts/functions/01_version_package_warnings.R"))

# load relevant functions
source(here("scripts/functions/03_get_taxon_id_function.R"))
source(here("scripts/functions/05_gbif_downstream_function.R"))
source(here("scripts/functions/06_itis_downstream_function.R"))

# add taxonomic information to the equation data.base

# define a function to get taxonomic information for each equation ID

# args
# equ.name.input - taxonomic name from the equation database
# equ.id - equation id from the equation database
# data_base - taxonomic database ( "itis" or "gbif" are supported)

get_taxonomic_info <- function(equ.name.input, equ.id, data_base = "itis") {
  
  # if the input name is a species then extract the genus
  z <- unlist( strsplit(x = equ.name.input, split = " ", fixed = TRUE) )
  if (length(z) > 1) {
    equ.name <- z[1]
  } else {equ.name <- equ.name.input}
  
  # get the taxon_id from the correct database
  taxon_id <- get_taxon_id(database_function = data_base, taxon_name = equ.name, ask_or_not = FALSE, tries = 5)
  
  # get the upwards classification for the taxon_id
  y.c <- classification(taxon_id)[[1]]
  y.c$ranknum <- 1:nrow(y.c)
  
  # get the taxonomic rank for the taxon_id
  equ.rank <- y.c[y.c$name == equ.name,][["ranknum"]]
  
  # get the taxonomic rank name for the taxon_id
  if (length(z) > 1) {
    equ.rank.name <- "species"
  } else { equ.rank.name <- y.c[y.c$name == equ.name,][["rank"]] }
  
  # get the order name and rank number for the taxon_id
  ord.rank <- y.c[y.c$rank == "order", ][["ranknum"]]
  ord.name <- y.c[y.c$rank == "order", ][["name"]]
  
  # if the taxon_id is above the rank of order, then the equation is unsuitable
  if( !("order" %in% y.c$rank) ) {
    
    message("equation above the rank of order")
    equ.suitable <- FALSE
    
  } else { equ.suitable <- TRUE }
  
  # stop the function if the equation is not suitable
  if (!equ.suitable) {
    
    # set the output for the unsuitable equation
    dmat.taxlist <- c(equation_id = equ.id, 
                      equation_name = equ.name.input, 
                      equation_rank = equ.rank.name,
                      tax_distance = NA,
                      tax_names = NA,
                      synonymns = NA)
    
    return(dmat.taxlist)
    
    stop("The equation is associated with a taxonomic level above order")
    
  }
  
  # get taxon id: database specific function i.e. get the taxon ID of the order of the taxa from the database
  ord.id <- get_taxon_id(database_function = data_base, 
                         taxon_name =  ord.name, 
                         ask_or_not = FALSE, tries = 5)
  ord.id <- ord.id[[1]]
  
  # get a distance matrix and list of taxa downstream of the order
  if (data_base == "gbif") {
    
    dmat.taxlist <- downstream_gbif(ord.id = ord.id, ord.name = ord.name)
    
  } else if (data_base == "itis") {
    
    dmat.taxlist <- downstream_itis(ord.id = ord.id, ord.name = ord.name)
    
    # add synonymns here
    equ.syn <- taxize::synonyms(taxon_id, db = data_base)
    equ.syn <- equ.syn[[1]]$syn_name
    if (is.null(equ.syn)) {
      equ.syn <- NA
    }
    
    # combine into dmat.taxlist
    dmat.taxlist <- c(dmat.taxlist, synonymns = equ.syn)
    
  } else {
    
    stop("Specify an appropriate database: itis or gbif")
    
  }
  
  # add the equation taxon name and the rank
  dmat.taxlist <- c(equation_id = equ.id, 
                    equation_name = equ.name.input, 
                    equation_rank = equ.rank.name,
                    dmat.taxlist)
  
  # return the list with relevant taxonomic information
  return(dmat.taxlist)
  
}

# test the function on a species
x <- get_taxonomic_info(equ.name.input = "Sinantherina socialis", 
                        equ.id = 31, 
                        data_base = "itis") 

# test the function on an order
y <- get_taxonomic_info(equ.name.input = "Insecta", 
                        equ.id = 45, 
                        data_base = "itis") 

# these functions give distance matrices to the genus level

# next step is to incorporate species level information
# to this...

# i.e. add one level of distance when it matches the genus

# this needs to be coupled to the equation database

# the search functions should be very simple

# but, will need some decisions e.g. if same distance but family then choose genus

# taxon name database
db.name <- "Pediciini"
db.name

# target name
tar.name <- "Limnophila costata"
tar.name

# this calculation across all the different distance matrices
d[which(row.names(d) == db.name), which(colnames(d) == tar.name) ]

# rank of all the different names



### END

