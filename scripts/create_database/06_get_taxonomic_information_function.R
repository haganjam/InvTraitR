
# Generate taxonomic information

# load relevant libraries
library(here)

# check for the correct packages:
source(here("scripts/create_database/01_version_package_warnings.R"))

# load relevant functions
source(here("scripts/functions/01_get_taxon_id_function.R"))
source(here("scripts/create_database/03_gbif_downstream_function.R"))
source(here("scripts/create_database/04_itis_downstream_function.R"))


# define a function to get the order for each taxon in the equation database

# args
# equ.name.input - taxonomic name from the equation database
# equ.id - equation id from the equation database
# data.base - taxonomic database ( "itis" or "gbif" are supported)
# life.stage - life stage of the equation

get_taxon_order <- function(equ.name.input, equ.id, data.base, life.stage = NA) {
  
  # if the input name is a species then extract the genus
  equ.name <- extract_genus(binomial = equ.name.input)
  
  # get the taxon_id from the correct database
  taxon_id <- get_taxon_id(database_function = data.base, taxon_name = equ.name, ask_or_not = FALSE, tries = 5)
  
  # if the taxon_id is not found in the database, then stop the function
  if (is.na(taxon_id)) {
    
    taxlist <- list(equation_id = equ.id, 
                 equation_name = equ.name.input,
                 equation_suitable = FALSE,
                 database = data.base,
                 equation_rank = NA,
                 order = NA,
                 life_stage = NA)
    
    return(taxlist)
    stop("Could not find this taxon name in the database")
    
  }
  
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
  
  # put this information into an output list
  taxlist <- list(equation_id = equ.id, 
                  equation_name = equ.name.input,
                  equation_suitable = equ.suitable,
                  database = data.base,
                  equation_rank = equ.rank.name,
                  order = ord.name,
                  life_stage = life.stage
                  )
  
  # add relevant synonymns
  if(data.base == "itis") {
    equ.syn <- taxize::synonyms(taxon_id, db = data.base)
    equ.syn <- equ.syn[[1]]$syn_name
    
    if (is.null(equ.syn)) {
      equ.syn <- NA
    }
    
  } else { equ.syn <- NA }
  
  # combine into dmat.taxlist
  taxlist <- c(taxlist, list(synonymns = equ.syn) )
  
  return(taxlist)
  
}

# define a function to get a taxonomic distance matrix from an order name

# args
# ord.name - name of the order
# data.base - taxonomic database ( "itis" or "gbif" are supported)

get_taxon_distance <- function(ord.name, data.base = "itis") {
  
  # get taxon id: database specific function i.e. get the taxon ID of the order of the taxa from the database
  ord.id <- get_taxon_id(database_function = data.base, 
                         taxon_name =  ord.name, 
                         ask_or_not = FALSE, tries = 5)
  ord.id <- ord.id[[1]]
  
  # if the taxon_id is not found in the database, then stop the function
  if (is.na(ord.id)) {
    
    dmat <- c(list(order = ord.name),
              list(tax_available = FALSE),
              list(tax_distance = NA,
                   tax_names = NA)) 
    
    return(dmat)
    stop("Could not find this order name in the database")
    
  }
  
  # get a distance matrix and list of taxa downstream of the order
  if (data.base == "gbif") {
    
    dmat <- downstream_gbif(ord.id = ord.id, ord.name = ord.name)
    
  } else if (data.base == "itis") {
    
    dmat <- downstream_itis(ord.id = ord.id, ord.name = ord.name)
    
  } else {
    
    stop("Specify an appropriate database: itis or gbif")
    
  }
  
  dmat <- c(list(order = ord.name),
            list(tax_available = TRUE),
            dmat)
  
  return(dmat)
  
}

### END
