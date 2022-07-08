#'
#' @title downstream_gbif()
#' 
#' @description Get downstream taxa from order from the gbif database
#' 
#' @details This function is used to get all downstream taxa from an order and process it
#' into a matrix of taxonomic distances and a vector of taxon names
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)
#' 
#' @param ord.id - taxon ID of the order
#' @param ord.name - name of the order
#' @param higher.tax.rank - rank of the ord.name (could be family or order)
#' 
#' @return list with two elements: 
#' 1. tax_distance - taxonomic distance matrix (saved as a sparse matrix); 
#' 2. tax_names - vector of all the names in the taxonomic distance matrix
#' 

downstream_gbif <- function(ord.id, ord.name, higher.tax.rank) {
  
  # load relevant packages
  library(dplyr)
  library(taxize)
  
  # get downstream with error handling and multiple tries
  source(here("scripts/create_database/02_get_downstream_taxa_function.R"))
  
  # if the higher taxonomic rank is order then
  if (higher.tax.rank == "order") {
    
    # get downstream taxa and process into a usable data.frame
    downtax.top <- get_downstream_taxa(sci_id = ord.id, downto = "family", db = "gbif", intermediate = FALSE)
    downtax.top <- downtax.top[[1]]
    downtax.top$parentname <- ord.name
    downtax.top$parentrank <- higher.tax.rank
    
    # loop over all families and get their downstream taxa
    downtax.int <- vector("list", length = nrow(downtax.top))
    for (i in 1:nrow(downtax.top)) {
      
      x <- get_downstream_taxa(sci_id = downtax.top$key[i], downto = "genus", db = "gbif", intermediate = FALSE)
      if (length(x[[1]]) == 0 ) {
        x <- NULL
      } else {
        x <- x[[1]]
        x$parentname <- downtax.top$name[i]
        x$parentrank <- "family"
        downtax.int[[i]] <- x
      }
    }
    
    # bind the intermediate taxa into a single data.frame
    downtax.top <- bind_rows(bind_rows(downtax.top), (downtax.int))
    
  } else if (higher.tax.rank == "family") {
    
    downtax.top <- get_downstream_taxa(sci_id = ord.id, downto = "genus", db = "gbif", intermediate = FALSE)
    downtax.top <- downtax.top[[1]]
    downtax.top$parentname <- higher.tax.name
    downtax.top$parentrank <- higher.tax.rank
    
  }
  
  # load the taxonomic distance matrix
  source(here("scripts/create_database/05_taxon_matrix_gbif_function.R"))
  
  # apply taxonomic weights
  weights <- mapply(function(x, y) tax.gbif[which(row.names(tax.gbif) == x), which(colnames(tax.gbif) == y) ],
                    x = downtax.top$rank,
                    downtax.top$parentrank)
  downtax.top$weights <- unlist(weights, use.names = FALSE)
  
  # create the distance matrix
  d.mat <- 
    downtax.top %>%
    select(from = parentname, to = name, weights)
  
  # use igraph to create a graph from the matrix
  d.g <- graph_from_data_frame(d = d.mat, directed=FALSE)
  
  # produce a distance matrix using the taxonomic weights
  d.g.dist <- distances(
    d.g,
    v = V(d.g),
    to = V(d.g),
    mode = c("all"),
    algorithm = c("bellman-ford")
  )
  
  # convert symmetrical values in upper matrix to zeros
  d.g.dist [upper.tri(d.g.dist , diag = FALSE)] <- 0
  
  # convert the distance matrix into a sparse matrix
  library(Matrix)
  d.g.dist  <- Matrix(d.g.dist, sparse = TRUE)
  
  # get a data.frame of unique taxonomic names and their 
  par.dat <- downtax.top[, c("parentname", "parentrank")]
  names(par.dat) <- c("name", "rank")
  tax.dat <- downtax.top[, c("name", "rank")]
  names(tax.dat) <- c("name", "rank")
  
  tax.names <- unique( rbind(par.dat, tax.dat) )
  
  # make a list of the distance matrix and the taxonomic names
  return(list(tax_distance = d.g.dist, tax_names = tax.names))
  
  }

### END
